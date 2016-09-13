///=== This is the Linear Regression for 4 helix parameters track fit algorithm.

///=== Written by: Thomas Schuh

#include "TMTrackTrigger/TMTrackFinder/interface/LinearRegression.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include <vector>
#include <algorithm>
#include <limits>

template< typename T >
std::ostream& operator << ( std::ostream& lhs, const std::pair< T, T >& rhs ) {
  return lhs << rhs.first << " " << rhs.second;
}

template< typename T >
float sgn( T val ) {
    return ( T( 0 ) < val ) - ( val < T( 0 ) );
}

std::pair< int, int > operator - ( std::pair< unsigned int, unsigned int >& lhs, std::pair< unsigned int, unsigned int >& rhs ) {
  return std::make_pair( std::abs( (int)lhs.first - (int)rhs.first ), std::abs( (int)lhs.second - (int)rhs.second ) );
}

template< typename T >
bool operator < ( std::pair< T, T >& lhs, std::pair< T, T >& rhs ) {
  if ( lhs.first < rhs.first and lhs.second < rhs.second )
    return true;
  return false;
}

void LinearRegression::initRun() {

  invPtToDphi_ = settings_->invPtToDphi();
  numPhiSectors_ = settings_->numPhiSectors();
  chosenRofPhi_ = settings_->chosenRofPhi();
  chosenRofZ_ = settings_->chosenRofZ();
  etaRegions_ = settings_->etaRegions();
  numFittingIterations_ = settings_->numTrackFitIterations();
  generalResidualCut_ = settings_->generalResidualCut();
  killingResidualCut_ = settings_->killingResidualCut();
  combineResiduals_ = settings_->combineResiduals();
  lineariseStubPosition_ = settings_->lineariseStubPosition();
  checkSectorConsistency_ = settings_->checkSectorConsistency();
  checkHTCellConsistency_ = settings_->checkHTCellConsistency();
  minStubLayers_ = settings_->minStubLayers();
  minPSLayers_ = settings_->minPSLayers();
  minPtToReduceLayers_ = settings_->minPtToReduceLayers();
  debug_ = settings_->debug() == 7;

};

L1fittedTrack LinearRegression::fit( const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg ) {

  initFit( l1track3D, iPhiSec, iEtaReg );
  if ( not checkValidity() )
    return createTrack( l1track3D );
  calcHelix( false );
  calcHelix();
  calcResidual();
  for ( int i = 1; i < numFittingIterations_ + 1; i++ ) {
    if ( not killLargestResidual() )
      break;
    updateLayerMap();
    if ( not checkValidity() )
      break;
    calcHelix();
    calcResidual();
  }
  return createTrack( l1track3D );

}

void LinearRegression::initFit( const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg ) {

  iPhiSec_ = iPhiSec;
  iEtaReg_ = iEtaReg;
  phiSectorSize_ = 2 * M_PI / float( numPhiSectors_ );
  phiCentre_ = phiSectorSize_ * ( 0.5 + float( iPhiSec_ ) ) - M_PI;
  zCentre_ = chosenRofZ_ * ( std::sinh( etaRegions_[ iEtaReg_ ] ) + std::sinh( etaRegions_[ iEtaReg_ + 1 ] ) ) / 2.; // based on asinh(eta) = cotan(theta)
  zSectorSize_ = chosenRofZ_ * std::fabs( std::sinh( etaRegions_[ iEtaReg_ ] ) - std::sinh( etaRegions_[ iEtaReg_ + 1 ] ) );
  houghNbinsPt_ = settings_->houghNbinsPt();
  houghNbinsPhi_ = settings_->houghNbinsPhi();
  binSizeQoverPtAxis_ = 2. / settings_->houghMinPt() / float( houghNbinsPt_ ) * invPtToDphi_;
  binSizePhiTrkAxis_ = 2. * M_PI / float( numPhiSectors_ ) / float( houghNbinsPhi_ );
  HTCell_ = l1track3D.getCellLocationRphi();
  if ( l1track3D.pt() > minPtToReduceLayers_ )
    minStubLayers_ --;
  stubs_ = l1track3D.getStubs();
  updateLayerMap();
  helixRPhi_ = l1track3D.getHelixRphi();
  helixRPhi_.first *= invPtToDphi_;
  helixRZ_ = std::make_pair( 0., zCentre_ / chosenRofZ_ );

}

void LinearRegression::updateLayerMap() {

  NumStubs_ = stubs_.size();
  nLayers_ = Utility::countLayers( settings_, stubs_ );
  nPSLayers_ = Utility::countLayers( settings_, stubs_, false, true );
  layerMap_.clear();
  psLayerMap_.clear();
  for( unsigned int is = 0; is < NumStubs_; is++ ) {
    const Stub* s( stubs_[ is ] );
    layerID_ = s->layerIdReduced();
    layerMap_[ layerID_ ].push_back( is );
    if ( s->psModule() )
      psLayerMap_[ layerID_ ].push_back( is );
  }

}

bool LinearRegression::checkValidity() {

  valid_ = true;
  if ( nLayers_ < minStubLayers_ )
    valid_ = false;
  if ( nPSLayers_ < minPSLayers_ )
    valid_ = false;
  return valid_;

}

void LinearRegression::calcHelix( bool useHelix ) {

  resetSums();
  for ( auto l : layerMap_ )  if ( useHelix or psLayerMap_[ l.first ].size() > 0 ) {
    phiMinMax_ = std::make_pair( std::numeric_limits< float >::infinity(), - std::numeric_limits< float >::infinity() );
    rMinMax_ = std::make_pair( std::numeric_limits< float >::infinity(), - std::numeric_limits< float >::infinity() );
    zMinMax_ = std::make_pair( std::numeric_limits< float >::infinity(), - std::numeric_limits< float >::infinity() );
    for ( auto is : l.second ) {
      const Stub* s( stubs_[ is ] );
      r_ = s->r();
      z_ = s->r() * helixRZ_.second + helixRZ_.first;
      if ( not s-> barrel() ) {
        r_ = ( s->z() - helixRZ_.first ) / helixRZ_.second;
        z_ = s->z();
      }
      phi_ = s->phi();
      if ( lineariseStubPosition_ )
        phi_ = reco::deltaPhi( s->phi(), reco::deltaPhi( r_ * helixRPhi_.first, std::asin( r_ * helixRPhi_.first ) ) );
      if ( not useHelix ) {
        r_ = s->r();
        phi_ = s->phi();
        z_ = s->z();
        if ( not s->psModule() )
          continue;
      }
      phi_ = reco::deltaPhi( phi_, phiCentre_ );
      z_ -= zCentre_;
      rMinMax_ = std::make_pair( std::min( rMinMax_.first, r_ ), std::max( rMinMax_.second, r_ ) );
      phiMinMax_ = std::make_pair( std::min( phiMinMax_.first, phi_ ), std::max( phiMinMax_.second, phi_ ) );
      zMinMax_ = std::make_pair( std::min( zMinMax_.first, z_ ), std::max( zMinMax_.second, z_ ) );
    }
    r_ = ( rMinMax_.first + rMinMax_.second ) / 2.;
    rTPhi_ = r_ - chosenRofPhi_;
    rTZ_ = r_ - chosenRofZ_;
    phi_ = ( phiMinMax_.first + phiMinMax_.second ) / 2.;
    z_ = ( zMinMax_.first + zMinMax_.second ) / 2.;
    updateSums();
  }
  calcLinearParameter();
  helixRPhi_ = std::make_pair( - slopeRPhi_, reco::deltaPhi( phiCentre_, -reco::deltaPhi( interceptRPhi_, slopeRPhi_ * chosenRofPhi_ ) ) );
  helixRZ_ = std::make_pair( zCentre_ + interceptRZ_ - slopeRZ_ * chosenRofZ_, slopeRZ_ );

}

void LinearRegression::calcResidual() {

  residRPhi_.clear();
  residRZ_.clear();
  for( unsigned int is = 0; is < NumStubs_; is++ ) {
    const Stub* s( stubs_[ is ] );
    r_ = s->r(); 
    phi_ = s->phi();
    if ( lineariseStubPosition_ )
      phi_ = reco::deltaPhi( s->phi(), reco::deltaPhi( r_ * helixRPhi_.first, std::asin( r_ * helixRPhi_.first ) ) );
    z_ = s->r() * helixRZ_.second + helixRZ_.first;
    if ( not s-> barrel() ) {
      r_ = ( s->z() - helixRZ_.first ) / helixRZ_.second;
      phi_ = reco::deltaPhi( helixRPhi_.second, r_ * helixRPhi_.first );
      z_ = s->z();
    }
    resid_ = reco::deltaPhi( reco::deltaPhi( helixRPhi_.second, helixRPhi_.first * r_ ), phi_ ) * r_ / s->sigmaPerp();
    if ( not s->barrel() )
      resid_ = reco::deltaPhi( s->phi(), phi_ ) * r_ / s->sigmaPerp();
    residRPhi_.push_back( resid_ );
    resid_ = ( s->z() - z_ ) / s->sigmaPar();
    if ( not s->barrel() )
      resid_ = ( s->r() - r_ ) / s->sigmaPar();
    residRZ_.push_back( resid_ );
  }

}

bool LinearRegression::killLargestResidual() {
 
  largestresid_ = -1.0;
  ilargestresid_  = -1;
  for ( auto l : layerMap_ ) if ( l.second.size() > 1 or NumStubs_ == nLayers_ ) for ( const unsigned int s : l.second ) {
    resid_ = std::max( residRPhi_[ s ], residRZ_[ s ] );
    if ( combineResiduals_ )
      resid_ = std::fabs( residRPhi_[ s ] ) / 2. + std::fabs( residRZ_[ s ] ) / 2.; 
    if ( resid_ > largestresid_ ) {
      largestresid_ = resid_;
      ilargestresid_ = s;
    }
  }
  if ( largestresid_ > killingResidualCut_ or NumStubs_ > nLayers_ ) {
    stubs_.erase( stubs_.begin() + ilargestresid_ );
    return true;
  }
  return false;

}

void LinearRegression::calcChiSq() {

  chiSq_ = 0.;
  for( unsigned int s = 0; s < NumStubs_; s++ )
    chiSq_ += residRPhi_[ s ] * residRPhi_[ s ] + residRZ_[ s ] * residRZ_[ s ];

}

L1fittedTrack LinearRegression::createTrack( const L1track3D& l1track3D ) {

  trackParams_["qOverPt"] = helixRPhi_.first / invPtToDphi_;
  trackParams_["phi0"] = helixRPhi_.second;
  trackParams_["z0"] =  helixRZ_.first;
  trackParams_["t"] = helixRZ_.second;
  if ( checkSectorConsistency_ ) {
    if ( std::fabs( helixRZ_.first - zCentre_ + helixRZ_.second * chosenRofZ_ ) > zSectorSize_ / 2. )
      valid_ = false;
    if ( std::fabs( reco::deltaPhi( reco::deltaPhi( helixRPhi_.second, phiCentre_ ), helixRPhi_.first * chosenRofPhi_ ) ) > phiSectorSize_ / 2. )
      valid_ = false;
    if ( std::fabs( trackParams_["qOverPt"] ) > 1. / settings_->houghMinPt() )
      valid_ = false;
  }
  if ( checkHTCellConsistency_ ) {
    helixRPhi_.second = reco::deltaPhi( reco::deltaPhi( helixRPhi_.second, phiCentre_ ), helixRPhi_.first * chosenRofPhi_ );
    trackCell_ = std::make_pair( (unsigned int)std::floor( (float)houghNbinsPt_ / 2. + helixRPhi_.first / binSizeQoverPtAxis_ ), (unsigned int)std::floor( (float)houghNbinsPhi_ / 2. + helixRPhi_.second / binSizePhiTrkAxis_ ) );
    if ( not ( ( trackCell_ - HTCell_ ) < std::make_pair( 3, 3 ) ) )
      valid_ = false;
  }
  return L1fittedTrack( settings_, l1track3D, stubs_, trackParams_["qOverPt"], 0, trackParams_["phi0"], trackParams_["z0"], trackParams_["t"], chiSq_, nPar_, iPhiSec_, iEtaReg_, valid_ );

}

void LinearRegression::resetSums() {

  N_ = 0;
  sumRPhi_ = 0.;
  sumRZ_ = 0.;
  sumRTPhi_ = 0.;
  sumRTZ_ = 0.;
  sumPhi_ = 0.;
  sumZ_ = 0.;
  sumRTPhi2_ = 0.;
  sumRTZ2_ = 0.;

}

void LinearRegression::updateSums() {

  N_ ++;
  sumRPhi_ += rTPhi_ * phi_;
  sumRZ_ += rTZ_ * z_;
  sumRTPhi_ += rTPhi_;
  sumRTZ_ += rTZ_;
  sumPhi_ += phi_;
  sumZ_ += z_;
  sumRTPhi2_ += rTPhi_ * rTPhi_;
  sumRTZ2_ += rTZ_ * rTZ_;

}

void LinearRegression::calcLinearParameter() {

  denominatorRPhi_ = N_ * sumRTPhi2_ - sumRTPhi_ * sumRTPhi_;
  denominatorRZ_ = N_ * sumRTZ2_ - sumRTZ_ * sumRTZ_;
  slopeRPhi_ = ( N_ * sumRPhi_ - sumRTPhi_ * sumPhi_ ) / denominatorRPhi_;
  slopeRZ_ = ( N_ * sumRZ_ - sumRTZ_ * sumZ_ ) / denominatorRZ_;
  interceptRPhi_ = ( sumRTPhi2_ * sumPhi_ - sumRTPhi_ * sumRPhi_ ) / denominatorRPhi_;
  interceptRZ_ = ( sumRTZ2_ * sumZ_ - sumRTZ_ * sumRZ_ ) / denominatorRZ_;

}