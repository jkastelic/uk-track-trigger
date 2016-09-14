#ifndef __HTCELL_H__
#define __HTCELL_H__

#include "TMTrackTrigger/TMTrackFinder/interface/Utility.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <vector>
#include <map>
#include <algorithm>
#include <utility>

class Settings;
class Stub;


//=== A single cell in a Hough Transform array.

class HTcell {

public:
  
  HTcell() {}
  ~HTcell() {}

  // Initialization with cfg params, boolean indicating if this is r-phi or r-z HT,   
  // rapidity range of current sector, estimated q/Pt of cell,
  // and (if called from r-phi HT) the bin number of the cell along the q/Pt axis of the r-phi HT array.
  void init(const Settings* settings, bool isRphiHT, float etaMinSector, float etaMaxSector, float qOverPt, unsigned int ibin_qOverPt = 0);

  // Add stub to this cell in HT array.
  void store (const Stub* stub) { vStubs_.push_back(stub); }

  // Add stub to this cell in HT array, indicating also which subsectors within the sector is consistent with.
  void store (const Stub* stub, const std::vector<bool>& inSubSecs) { this->store(stub); subSectors_[stub] = inSubSecs; if (inSubSecs.size() != numSubSecs_) throw cms::Exception("HTcell: Wrong number of subsectors!");}

  // Termination. Search for track in this HT cell etc.
  void end();

  //=== Get results
  //=== Most of these functions operate on the filtered stubs, which are stubs passing any requested stub filters (e.g. bend filter). 
  //=== If no filters were requested, they are identical to the unfiltered stubs.)

  // Get filtered stubs in this cell in HT array.
  const std::vector<const Stub*>& stubs() const { return vFilteredStubs_; }

  // Check if a specific stub is in this cell and survived filtering.
  bool stubInCell( const Stub* stub ) const { return (std::count(vFilteredStubs_.begin(), vFilteredStubs_.end(), stub ) > 0); }

  // Check if a specific stub was stored to this cell (without checking if it survived filtering).
  bool stubStoredInCell( const Stub* stub ) const { return (std::count(vStubs_.begin(), vStubs_.end(), stub ) > 0); }

  // Return info useful for deciding if there is a track candidate in this cell.
  unsigned int numStubs()         const { return vFilteredStubs_.size(); }      // Number of filtered stubs 
  unsigned int numLayers()        const { return numFilteredLayersInCell_; }    // Number of tracker layers with filtered stubs
  unsigned int numLayersSubSec()  const { return numFilteredLayersInCellBestSubSec_; }  // Number of tracker layers with filtered stubs,  requiring all stubs to be in same subsector to be counted. The number returned is the highest layer count found in any of the subsectors in this sector. If subsectors are not used, it is equal to numLayers().

  // Useful for debugging.
  unsigned int numUnfilteredStubs()   const { return vStubs_.size(); }    // Number of unfiltered stubs 

  //=== Check if stubs in this cell form valid track candidate.

  // N.B. If subsectors within a sector are not being used, then numFilteredLayersInCellBestSubSec_ = numFilteredLayersInCell_.
  // WARNING: If some tracks are killed as the r-phi HT array can't read them out within the TM period, 
  // killed tracks are still found by this function. It is in HTbase::calcTrackCands2D() that they are killed.
  bool trackCandFound() const { 
    return ( (fabs(qOverPtCell_) > 1/minPtToReduceLayers_)  ?  
             (numFilteredLayersInCellBestSubSec_ >= minStubLayers_)  :  (numFilteredLayersInCellBestSubSec_ >= minStubLayers_ - 1) );
  }

  //=== Disable filters (used for debugging).

  void disableBendFilter() {useBendFilter_ = false;}

private:
  // Calculate how many tracker layers the filtered stubs in this cell are in
  unsigned int calcNumFilteredLayers()  const { return Utility::countLayers( settings_, vFilteredStubs_ ); }

  // Calculate how many tracker layers the filter stubs in this cell are in, when only the subset of those stubs
  // that are in the specified subsector are counted.
  unsigned int calcNumFilteredLayers(unsigned int iSubSec) const;

  // Estimate track bend angle at a given radius, derived using the track q/Pt at the centre of this HT cell, ignoring scattering.
  float dphi(float rad) const { return (invPtToDphi_ * rad * qOverPtCell_); }

  // Produce a filtered collection of stubs in this cell that all have consistent bend
  std::vector<const Stub*> bendFilter( const std::vector<const Stub*>& stubs ) const;

  // Filter stubs so as to prevent more than specified number of stubs being stored in one cell.
  // This reflects finite memory of hardware.
  std::vector<const Stub*> maxStubCountFilter( const std::vector<const Stub*>& stubs ) const;

private:
  //=== Configuration parameters

  const Settings* settings_;

  // Note if if this r-phi or r-z HT.
  bool isRphiHT_;

  float etaMinSector_; // rapidity range of this sector.
  float etaMaxSector_;

  float qOverPtCell_; // track q/Pt corresponding to centre of this cell.

  // Note bin number of cell along q/Pt axis of r-phi HT array. (Not used for r-z HT).
  unsigned int ibin_qOverPt_;

  float invPtToDphi_; // B*c/2E11

  // Min. numbers of layers with stubs required to produce track candidate.
  float  minStubLayers_;
  // Reduce MinStubLayers by 1 for HT cells corresponding to track pt above this cut.
  // If this is set to >= 10000., this option is disabled.
  float  minPtToReduceLayers_;
 
  // Use filter in each HT cell using only stubs which have consistent bend?
  bool   useBendFilter_;
  // A filter is used each HT cell, which prevents more than the specified number of stubs being stored in the cell. (Reflecting memory limit of hardware).   
  unsigned int maxStubsInCell_;

  // Note if daisy-chain firmware is in use, together with digitized stubs.
  bool daisyChainFirmware_;

  // Number of subsectors (if any) within each sector.
  unsigned int numSubSecs_;

  //=== data

  std::vector<const Stub*> vStubs_; // Stubs in this cell
  std::vector<const Stub*> vFilteredStubs_; // Stubs in cell selected by applying all requested stub filters (e.g. bend and/or eta filter ...)

  unsigned int numFilteredLayersInCell_; // How many tracker layers these filtered stubs are in
  unsigned int numFilteredLayersInCellBestSubSec_; // Ditto, but requiring all stubs to be in same subsector to be counted. This number is the highest layer count found in any of the subsectors in this sector.

  std::map<const Stub*, std::vector<bool>> subSectors_; // Indicate which subsectors within the sector this stub is consistent with.
};
#endif

