#include "TMTrackTrigger/TMTrackFinder/interface/HTrz.h"
#include "TMTrackTrigger/TMTrackFinder/interface/InputData.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"

#include <vector>
#include <set>

//=== The r-z Hough Transform array for a single (eta,phi) sector.
//===
//=== Its axes are (z0, zTrk), where zTrk is the z at which the track crosses a 
//=== user-configurable radius from the beam-line.

using namespace std;

// Maximum |gradient| of line corresponding to any stub. Should be less than the value of 1.0 assumed by the firmware.
float        HTrz::maxLineGradient_ = 0.;
// Error count when stub added to cell which does not lie NE, E or SE of stub added to previous HT column.
unsigned int HTrz::numErrorsTypeA_ = 0;
// Error count when stub added to more than 2 cells in one HT column (problem only for Thomas' firmware).
unsigned int HTrz::numErrorsTypeB_ = 0;
// Error count normalisation
unsigned int HTrz::numErrorsNormalisation_ = 0;

//=== Initialise
 
void HTrz::init(const Settings* settings, float etaMinSector, float etaMaxSector, float qOverPt) {
  HTbase::settings_  = settings;

  int nCellsHT = settings->houghNcellsRz(); // Total number of required cells in HT array (if > 0)

  //--- Specification of HT z0 axis 

  maxAbsZ0Axis_  = settings->beamWindowZ(); // Half-width of z0 axis.
  nBinsZ0Axis_   = settings->houghNbinsZ0(); // no. of bins in HT array in z0.
  if (nCellsHT > 0) nBinsZ0Axis_ = 1; // Will calculate number of bins automatically. Initialize it to non-zero value.
  binSizeZ0Axis_ = 2*maxAbsZ0Axis_ / nBinsZ0Axis_;

  //--- Specification of HT zTrk axis. 

  // N.B. zTrk corresponds to the z coordinate where the track crosses radius = chosenRofZ.
  chosenRofZ_ = settings->chosenRofZ(); 
  // Although we talk of eta regions, actually these regions are defined according to the z coordinate
  // where the track intercepts radius = chosenRofZ_. (See Sector::init()).
  // This z coordinate is named zTrk.
  minZtrkAxis_ = chosenRofZ_ / tan( 2. * atan(exp(-etaMinSector)) );
  maxZtrkAxis_ = chosenRofZ_ / tan( 2. * atan(exp(-etaMaxSector)) );
  nBinsZtrkAxis_   = settings->houghNbinsZ65(); // no. of bins in HT array in z65 (or any other z related variable)
  if (nCellsHT > 0) nBinsZtrkAxis_ = 1; // Will calculate number of bins automatically. Initialize it to non-zero value.
  binSizeZtrkAxis_ = (maxZtrkAxis_ - minZtrkAxis_) / nBinsZtrkAxis_;

  // Did user specify number of cells required in HT array? If so, determine number of bins along
  // array axes such that their product equals required number of cells, and that their ratio gives
  // a maximum line |gradient| of stubs crossing the array of 1.0.
  if (nCellsHT > 0) { 
    // Get line gradient with current array axes.
    float currentLineGrad = this->calcMaxLineGradArray();
    // Calculate new number of bins on each axis to meet constraint.
    float fact = nBinsZ0Axis_ * currentLineGrad / nBinsZtrkAxis_;
    nBinsZ0Axis_   = ceil( sqrt(nCellsHT * fact) );
    nBinsZtrkAxis_ = int ( sqrt(nCellsHT / fact) );
    // And recalculate bin size accordingly.
    binSizeZ0Axis_   = 2*maxAbsZ0Axis_ / nBinsZ0Axis_;
    binSizeZtrkAxis_ = (maxZtrkAxis_ - minZtrkAxis_) / nBinsZtrkAxis_;
  }

  // Note max. |gradient| that the line corresponding to any stub in any of the r-phi HT arrays could have.
  // Firmware assumes this should not exceed 1.0;
  HTrz::maxLineGradient_ = max( HTrz::maxLineGradient_, this->calcMaxLineGradArray());

  //--- Other options used when filling the HT.

  // Don't fill all HT cells nominally crossed by line corresponding to stub.
  killSomeHTCellsRz_ = settings->killSomeHTCellsRz();
  // Should algorithm allow for uncertainty in stub (r,z) coordinate caused by length of 2S module strips when fill stubs in r-z HT?
  handleStripsRzHT_  = settings->handleStripsRzHT();

  //--- Options for duplicate track removal after running HT.
  unsigned int dupTrkAlgRz = settings->dupTrkAlgRz();
  HTbase::killDupTrks_.init(settings, dupTrkAlgRz);

  // Resize HT array to suit these specifications, and initialise each cell with configuration parameters.
  HTbase::htArray_.resize(nBinsZ0Axis_, nBinsZtrkAxis_, false);

  const bool isRphiHT = false;
  for (unsigned int i = 0; i < nBinsZ0Axis_; i++) {
    for (unsigned int j = 0; j < nBinsZtrkAxis_; j++) {
      HTbase::htArray_(i,j).init( settings, isRphiHT, etaMinSector, etaMaxSector, qOverPt ); // Calls HTcell::init()
    }
  }

  static set<float> first;
  if (std::count(first.begin(), first.end(), etaMinSector) == 0) {
    first.insert(etaMinSector);
    cout<<"=== R-Z HOUGH TRANSFORM AXES RANGES: abs(z0) < "<<maxAbsZ0Axis_<<" & "<<minZtrkAxis_<<" < zTrk < "<<maxZtrkAxis_<<" ==="<<endl<<endl;
    cout<<"=== R-Z HOUGH TRANSFORM ARRAY SIZE: z0 bins = "<<nBinsZ0Axis_<<" zTrk bins = "<<nBinsZtrkAxis_<<endl; 
  }
}

//=== Add stub to HT array.

void HTrz::store( const Stub* stub) {


  // Loop over z0 related bins in HT array.
  for (unsigned int i = 0; i < nBinsZ0Axis_; i++) {

    // In this tan_lambda bin, find the range of zTrk bins that this stub is consistent with.
    pair<unsigned int, unsigned int> iRange = this->iZtrkRange( stub, i);
    unsigned int iZtrkBinMin = iRange.first;
    unsigned int iZtrkBinMax = iRange.second;

    // Store stubs in these cells.
    for (unsigned int j = iZtrkBinMin; j <= iZtrkBinMax; j++) {  
      HTbase::htArray_(i, j).store( stub ); // Calls HTcell::store()
    }

    // Check that limitations of firmware would not prevent stub being stored correctly in this HT column.
    this->countFirmwareErrors(i, iZtrkBinMin, iZtrkBinMax);
  }
}

//=== For a given z0 bin, find the range of zTrk bins that a given stub is consistent with.
//=== Return as a pair (min bin, max bin)
//=== If it range lies outside the HT array, then the min bin will be set larger than the max bin.

pair<unsigned int, unsigned int> HTrz::iZtrkRange( const Stub* stub, unsigned int iZ0Bin, bool debug) const {

  // Note z0 value corresponding to centre of this bin.
  float z0Bin    = -maxAbsZ0Axis_ + (iZ0Bin + 0.5) * binSizeZ0Axis_;
  // Note change in this z0 value needed to reach either edge of the bin. 
  float z0BinVar = 0.5*binSizeZ0Axis_;

  // Calculate range of zTrk that would allow a track in this z0 range to pass through the stub.
  float zTrk    = ( stub->z() * chosenRofZ_ + z0Bin    *     (stub->r() - chosenRofZ_) ) / stub->r();
  float zTrkVar = (                           z0BinVar * fabs(stub->r() - chosenRofZ_) ) / stub->r();
  float zTrkMin = zTrk - zTrkVar;
  float zTrkMax = zTrk + zTrkVar;

  // Allow for uncertainty due to strip length if requested.
  if (handleStripsRzHT_) {
    // Estimate uncertainty due to strip length, using first order derivative of zTrk w.r.t. stub coords.
    // Note that barrel modules only care about zErr and endcap ones about rErr.
    float zTrkVarStub;
    if (stub->barrel()) {
      zTrkVarStub = stub->zErr() * chosenRofZ_ / stub->r();
    } else {
      zTrkVarStub = fabs(stub->z() - z0Bin) * (chosenRofZ_ / stub->r()) * (stub->rErr() / stub->r() );
    }
    //if (zTrkVarStub > 7) cout<<"CHECKRES "<<zTrkVarStub<<" "<<zTrkVar<<" "<<stub->barrel()<<" "<<stub->r()<<" "<<stub->z()<<" "<<stub->eta()<<endl;
    zTrkMin -= zTrkVarStub; 
    zTrkMax += zTrkVarStub; 
  }

  // Allow for multiple scattering/resolution
  // zTrkMin -= 0.005;
  // zTrkMax += 0.005;

  pair<float, float> zTrkRange( zTrkMin, zTrkMax );

  // Determine which HT array cell range in track-phi this range "z0Range" corresponds to.
  pair<unsigned int, unsigned int> iZtrkBinRange = this->HTbase::convertCoordRangeToBinRange(zTrkRange, nBinsZtrkAxis_, minZtrkAxis_, binSizeZtrkAxis_, killSomeHTCellsRz_);

  return iZtrkBinRange;
}

//=== Check that limitations of firmware would not prevent stub being stored correctly in this HT column.

void HTrz::countFirmwareErrors(unsigned int iZ0Bin, unsigned int iZtrkBinMin, unsigned int iZtrkBinMax) {
  static unsigned int iZtrkBinMinLast = 0;
  static unsigned int iZtrkBinMaxLast = 99999;
  // Reinitialize if this is left-most column in HT array.
  if (iZ0Bin == 0) {
    iZtrkBinMinLast = 0;
    iZtrkBinMaxLast = 99999;
  }

  // Only do check if stub is being stored somewhere in this HT column.
  if (iZtrkBinMax >= iZtrkBinMin) {
    //--- Remaining code below checks that firmware could successfully store this stub in this column.
    //   (a) Does cell lie NE, E or SE of cell filled in previous column?
    bool OK_a = (iZtrkBinMin + 1 >= iZtrkBinMinLast) && (iZtrkBinMax <= iZtrkBinMaxLast + 1);
    //   (b) Are no more than 2 cells filled in this column (problem only for Thomas' firmware)
    bool OK_b = (iZtrkBinMax - iZtrkBinMin + 1 <= 2);

    if ( ! OK_a ) numErrorsTypeA_++;
    if ( ! OK_b ) numErrorsTypeB_++;
    numErrorsNormalisation_++; // No. of times a stub is added to an HT column.

    iZtrkBinMinLast = iZtrkBinMin;
    iZtrkBinMaxLast = iZtrkBinMax;
  }
}

//=== Get the values of the track helix params corresponding to middle of a specified HT cell (i,j).
//=== The helix parameters returned will be those corresponding to the two axes of the HT array.
//=== So they might be (z0, zTrk), (z0, tan_lambda), (z0, eta) etc. depending on the configuration.

pair<float, float> HTrz::helix2Dhough(unsigned int i, unsigned int j) const {
  float z0   = -maxAbsZ0Axis_ + (i + 0.5) * binSizeZ0Axis_;
  float zTrk = minZtrkAxis_   + (j + 0.5) * binSizeZtrkAxis_;
  return pair<float, float>(z0, zTrk); 
}

//=== Get the values of the track helix params corresponding to middle of a specified HT cell (i,j).
//=== The helix parameters returned will be always be (z0, tan_lambda), irrespective of how the axes
//=== of the HT array are defined.

pair<float, float> HTrz::helix2Dconventional(unsigned int i, unsigned int j) const {
  // Get the helix parameters corresponding to the axes definitions of the HT.
  pair<float, float> helix2Dht = this->helix2Dhough(i,j);
  // Convert to the conventionally agreed pair of helix parameters, (z0, tan_lambda).
  float z0   = helix2Dht.first; // easy
  float zTrk = helix2Dht.second; 
  float tanLambda = (zTrk - z0)/chosenRofZ_;
  return pair<float, float>(z0, tanLambda);
}

//=== Which cell in HT array should this TP be in, based on its true trajectory?
//=== Returns (-1,-1) if TP not expected to be in any cell in this array.

pair<int, int> HTrz::trueCell( const TP* tp ) const {

  // Get HT axis variables corresponding to this TP.
  float z0    = tp->z0();
  float zTrk  = tp->trkZAtR( chosenRofZ_ );
  // Convert to bin numbers inside HT array.
  int iZ0   = floor( (  z0   - ( -maxAbsZ0Axis_) ) / binSizeZ0Axis_   );
  int iZtrk = floor( (  zTrk -    minZtrkAxis_   ) / binSizeZtrkAxis_ );
  if (iZ0 >= 0 && iZ0 < int(nBinsZ0Axis_) && iZtrk >= 0 && iZtrk < int(nBinsZtrkAxis_)) {
    return pair<int, int>(iZ0, iZtrk); // Cell found, so return it.
  } else {
    return pair<int, int>(-1, -1); // TP is not in this HT array at all.
  }
}

//=== Which cell in HT array should this fitted track be in, based on its fitted trajectory?
//=== Returns (-1,-1) if fitted track not expected to be in any cell in this array.

pair<int, int> HTrz::getCell( const L1fittedTrack* fitTrk ) const {

  // Get HT axis variables corresponding to this TP.
  float z0    = fitTrk->z0();
  // Convert (z0, tanLambda) of track at chosen radius used by HT.
  float zTrk  = z0 + chosenRofZ_ * fitTrk->tanLambda(); // neglects transverse impact parameter & track curvature.
  // Convert to bin numbers inside HT array.
  int iZ0   = floor( (  z0   - ( -maxAbsZ0Axis_) ) / binSizeZ0Axis_   );
  int iZtrk = floor( (  zTrk -    minZtrkAxis_   ) / binSizeZtrkAxis_ );
  if (iZ0 >= 0 && iZ0 < int(nBinsZ0Axis_) && iZtrk >= 0 && iZtrk < int(nBinsZtrkAxis_)) {
    return pair<int, int>(iZ0, iZtrk); // Cell found, so return it.
  } else {
    return pair<int, int>(-1, -1); // TP is not in this HT array at all.
  }
}

//=== Calculate maximum |gradient| that any stub's line across this HT array could have, so can check it doesn't exceed 1.

float HTrz::calcMaxLineGradArray() const {
  // Get max. |gradient| possible in this HT array.
  // N.B. This ignores the effect of option handleStripsRzHT_, which may not matter?
  float gradOuter = fabs(1.0 - chosenRofZ_/settings_->trackerOuterRadius());
  float gradInner = fabs(1.0 - chosenRofZ_/settings_->trackerInnerRadius());
  float maxGrad = max(gradOuter, gradInner);
  // Convert it to units of bin width.
  maxGrad *= binSizeZ0Axis_/binSizeZtrkAxis_;
  return maxGrad;
}

