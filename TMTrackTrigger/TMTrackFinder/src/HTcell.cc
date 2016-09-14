#include <TMTrackTrigger/TMTrackFinder/interface/HTcell.h>
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"

//=== Initialization with cfg params, boolean indicating if this is r-phi or r-z HT, 
//=== rapidity range of current sector, and estimated q/Pt of cell,
//=== and (if called from r-phi HT) the bin number of the cell along the q/Pt axis of the r-phi HT array.

void HTcell::init(const Settings* settings, bool isRphiHT, float etaMinSector, float etaMaxSector, float qOverPt, unsigned int ibin_qOverPt) {

  settings_ = settings;

  // Note if if this r-phi or r-z HT.
  isRphiHT_ = isRphiHT;

  // Note track q/Pt. 
  // In this case of an r-phi HT, each cell corresponds to a unique q/Pt.
  // In the case of an r-z HT, it is assumed that we know q/Pt from previously run r-phi HT.
  qOverPtCell_ = qOverPt;
  // Note bin number of cell along q/Pt axis of r-phi HT array. (Not used if r-z HT).
  ibin_qOverPt_ = ibin_qOverPt;
  // Rapidity range of sector.
  etaMinSector_ = etaMinSector;
  etaMaxSector_ = etaMaxSector;

  // Min. number of layers in HT cell that must have stubs for track to be declared found.
  minStubLayers_ = settings->minStubLayers();
  // Reduce MinStubLayers by 1 for HT cells corresponding to track pt above this cut.
  // If this is set to >= 10000., this option is disabled.
  minPtToReduceLayers_ = settings->minPtToReduceLayers();

  // The following only relevant to r-phi Hough transform.
  if (isRphiHT_) {
    invPtToDphi_   = settings->invPtToDphi();  // B*c/2E11

    // Use filter in each HT cell using only stubs which have consistent bend?
    useBendFilter_ = settings->useBendFilter();
  }

  // A filter is used each HT cell, which prevents more than the specified number of stubs being stored in the cell. (Reflecting memory limit of hardware).
  maxStubsInCell_ = settings->maxStubsInCell();

  // Note if daisy-chain firmware is in use, together with digitized stubs.
  daisyChainFirmware_ = (settings->firmwareType() == 1) && settings->enableDigitize();

  // Check if subsectors are being used within each sector. These are only ever used for r-phi HT.
  numSubSecs_ = isRphiHT_   ?   settings->numSubSecsEta()  :  1;
}

//=== Termination. Search for track in this HT cell etc.

void HTcell::end(){
  // Produce list of filtered stubs by applying all requested filters (e.g. on stub bend).
  // (If no filters are requested, then filtered & unfiltered stub collections will be identical).

  // N.B. Other filters,  such as the r-z filters, which the firmware runs after the HT because they are too slow within it,
  // are not defined here, but instead inside class TrkFilterAfterRphiHT.

  vFilteredStubs_ = vStubs_;
  // The bend filter is only relevant to r-phi Hough transform.
  if (isRphiHT_) {
    if (useBendFilter_) vFilteredStubs_ = this->bendFilter(vFilteredStubs_);
  }
  // Prevent too many stubs being stored in a single HT cell if requested (to reflect hardware memory limits).
  // N.B. This MUST be the last filter applied.
  if (maxStubsInCell_ <= 99) vFilteredStubs_ = this->maxStubCountFilter(vFilteredStubs_);

  // Calculate the number of layers the filtered stubs in this cell are in.
  numFilteredLayersInCell_ = this->calcNumFilteredLayers();

  if (numSubSecs_ > 1) { 
    // If using subsectors within each sector, calculate the number of layers the filters stubs in this cell are in,
    // when one considers only the subset of the stubs within each subsector.
    // Look for the "best" subsector.
    numFilteredLayersInCellBestSubSec_ = 0;
    for (unsigned int i = 0; i < numSubSecs_; i++) {
      unsigned int numLaySubSec = this->calcNumFilteredLayers(i);
      numFilteredLayersInCellBestSubSec_ = std::max(numFilteredLayersInCellBestSubSec_, numLaySubSec);
    }
  } else {
    // If only 1 sub-sector, then subsector and sector are identical.
    numFilteredLayersInCellBestSubSec_ = numFilteredLayersInCell_;
  }
}

// Calculate how many tracker layers the filter stubs in this cell are in, when only the subset of those stubs
// that are in the specified subsector are counted.

unsigned int HTcell::calcNumFilteredLayers(unsigned int iSubSec) const {
  std::vector<const Stub*> stubsInSubSec;
  for (const Stub* s : vFilteredStubs_) {
    const std::vector<bool>& inSubSec = subSectors_.at(s); // Find out which subsectors this stub is in.
    if (inSubSec[iSubSec]) stubsInSubSec.push_back(s);
  }
  return Utility::countLayers( settings_, stubsInSubSec );
}


//=== Produce a filtered collection of stubs in this cell that all have consistent bend.
//=== Only called for r-phi Hough transform.

std::vector<const Stub*> HTcell::bendFilter( const std::vector<const Stub*>& stubs ) const
{
	using namespace std;
	
  // Create bend-filtered stub collection.
  vector<const Stub*> filteredStubs;
  for (const Stub* s : stubs) {

    // Require stub bend to be consistent with q/Pt of this cell.

    if (daisyChainFirmware_) {
      // Daisy chain firmware doesn't have access to variables needed to calculate dphi of stub,
      // but instead knows integer range of q/Pt bins that stub bend is compatible with, so use these.
      if (s->min_qOverPt_bin() <= ibin_qOverPt_ && ibin_qOverPt_ <= s->max_qOverPt_bin() )  filteredStubs.push_back(s);
    } else {
      // Systolic array & 2-c-bin firmware do hace access to stub dphi, so can use it.
      // Predict track bend angle based on q/Pt of this HT cell and radius of stub.
      float predictedDphi = this->dphi( s->r() );
      // Require reconstructed and predicted values of this quantity to be consistent within estimated resolution. 
      if (fabs(s->dphi() - predictedDphi) < s->dphiRes()) filteredStubs.push_back(s);
    }
  }
  return filteredStubs;
}

//=== Filter stubs so as to prevent more than specified number of stubs being stored in one cell.
//=== This reflects finite memory of hardware.

std::vector<const Stub*> HTcell::maxStubCountFilter( const std::vector<const Stub*>& stubs ) const
{
	using namespace std;
	
  vector<const Stub*> filteredStubs;
  unsigned int numStubsToDelete = (stubs.size() > maxStubsInCell_)  ?  stubs.size() - maxStubsInCell_  :  0;
  // If there are too many stubs in a cell, the hardware throws away the first ones and keeps the last ones. 
  for (unsigned int i = numStubsToDelete; i < stubs.size(); i++) {
    filteredStubs.push_back(stubs[i]);
  }
  return filteredStubs;
}
