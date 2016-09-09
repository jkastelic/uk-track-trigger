///=== This is the Linear Regression for 4 helix parameters track fit algorithm.
 
///=== Written by: Thomas Schuh
 
#ifndef __LINEARREGRESSION__
#define __LINEARREGRESSION__

#include "TMTrackTrigger/TMTrackFinder/interface/TrackFitGeneric.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"

#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include <vector>
#include <sstream>
#include <string>

class LinearRegression : public TrackFitGeneric {

public:

    LinearRegression(const Settings* settings, const uint nPar) : TrackFitGeneric( settings ), settings_( settings ), nPar_( nPar ) {};
 
    virtual ~LinearRegression() {};
 
    virtual void initRun();

    L1fittedTrack fit(const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg);
 
    std::string getParams() { return "LinearRegression"; };
 
protected:

  void initFit( const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg );
  void updateLayerMap();
  bool checkValidity();
  void calcHelix( bool useHelix = true );
  void calcResidual();
  bool killLargestResidual();
  void calcChiSq();
  L1fittedTrack createTrack( const L1track3D& l1track3D );
  void resetSums();
  void updateSums();
  void calcLinearParameter();

  // LinearRegression(const Settings* settings, const uint nPar)
  const Settings* settings_;
  uint nPar_;

  // initRun()
  float invPtToDphi_;
  unsigned int numPhiSectors_;
  std::vector< double > etaRegions_;
  float chosenRofPhi_;
  float chosenRofZ_;
  int numFittingIterations_;
  double generalResidualCut_;
  double killingResidualCut_;
  bool combineResiduals_;
  bool lineariseStubPosition_;
  bool checkSectorConsistency_;
  bool checkHTCellConsistency_;
  unsigned int minPSLayers_;
  unsigned int minStubLayers_;
  float minPtToReduceLayers_;
  bool debug_;

  // initFit( const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg )
  std::stringstream debugStream_;
  unsigned int iPhiSec_;
  unsigned int iEtaReg_;
  float phiCentre_;
  float zCentre_;
  float zSectorSize_;
  float phiSectorSize_;
  std::vector<const Stub*> stubs_;
  std::pair< float, float > helixRPhi_;
  std::pair< float, float > helixRZ_;
  unsigned int houghNbinsPt_;
  unsigned int houghNbinsPhi_;
  float binSizeQoverPtAxis_;
  float binSizePhiTrkAxis_;
  std::pair< unsigned int, unsigned int > HTCell_;
  std::pair< unsigned int, unsigned int > trackCell_;

  // updateLayerMap()
  unsigned int NumStubs_;
  unsigned int layerID_;
  std::map< unsigned int, std::vector< unsigned int > > layerMap_;
  std::map< unsigned int, std::vector< unsigned int > > psLayerMap_;
  unsigned int nLayers_;
  unsigned int nPSLayers_;

  // checkValidity( const L1track3D& l1track3D )
  bool valid_;

  // calcHelixRPhi(), calcHelixRZ()
  std::pair< float, float > rMinMax_;
  std::pair< float, float > phiMinMax_;
  std::pair< float, float > zMinMax_;
  float r_;
  float rTPhi_;
  float rTZ_;
  float phi_;
  float z_;
  float slopeRPhi_;
  float slopeRZ_;
  float interceptRPhi_;
  float interceptRZ_;
 
  // calcResidualRPhi(), calcResidualRZ()
  float rHelix_;
  float phiHelix_;
  float moduleRadian_;
  float modulePhi_;
  float resid_;
  std::vector< float > residRPhi_;
  std::vector< float > residRZ_;

  // killLargestResidual()
  float largestresid_;
  int ilargestresid_;
 
  // calcChiSq()
  float chiSq_;

  // createTrack( const L1track3D& l1track3D )
  std::map< std::string, double > trackParams_;

  // resetSums(), updateSums()
  unsigned int N_;
  float sumRPhi_, sumRZ_, sumRTPhi_, sumRTZ_, sumPhi_, sumZ_, sumRTPhi2_, sumRTZ2_;

  // calcLinearParameter()
  float denominatorRPhi_;
  float denominatorRZ_;

};
#endif