#ifndef __TrkRZfilter_H__
#define __TrkRZfilter_H__

#include "TMTrackTrigger/TMTrackFinder/interface/L1track2D.h"
#include "TMTrackTrigger/TMTrackFinder/interface/HTrphi.h"

#include <vector>

class Settings;
class Stub;



//=== This class runs filters in track candidates previously found by the r-phi Hough transform
//=== e.g. Filters requiring the stubs to be consistent with a straight line in the r-z plane.
//===
//=== The filtering removes inconsistent stubs from the track candidates, & also kills some track candidates
//=== altogether if the filter leaves them with too few stubs.
//===
//=== Some r-z filters may also add an estimate of the r-z helix parameters to the selected track candidates.
//===
//=== It does NOT contain filters such as the bend filter, which are so simple that the firmware can run them 
//=== INSIDE the r-phi HT. Simple filters of this kind are in class HTcell.

class TrkRZfilter {

public:

  TrkRZfilter() {}
  ~TrkRZfilter() {}
  
  // Initialize configuration parameters, and note eta range covered by sector and phi coordinate of its centre.
  void init(const Settings* settings, float etaMinSector, float etaMaxSector, float phiCentreSector);

  // Filters track candidates (found by the r-phi Hough transform), removing inconsistent stubs from the tracks, 
  // also killing some of the tracks altogether if they are left with too few stubs.
  // Also adds an estimate of r-z helix parameters to the selected track objects, if the filters used provide this.
  std::vector<L1track2D> filterTracks(const std::vector<L1track2D>& tracks);

  //=== Extra information about each track input to filter. (Only use after you have first called filterTracks).

  // Number of seed combinations considered by the ZTrk Filter for each input track.
  std::vector<unsigned int> numZtrkSeedCombsPerTrk() const {return numZtrkSeedCombsPerTrk_; }
  // Number of seed combinations considered by the Seed Filter for each input track.
  std::vector<unsigned int> numSeedCombsPerTrk() const {return numSeedCombsPerTrk_; }
  std::vector<unsigned int> numGoodSeedCombsPerTrk() const {return numGoodSeedCombsPerTrk_; } // Only counts seeds compatible with beam-spot.


private:

  // Check if a track candidate with stubs in specified number of tracker layers & given estimated q/Pt has stubs in enough
  // layers to be defined as a valid track candidate.
  bool trackCandCheck(unsigned int nLayers, float trkQoverPt) const { 
    return ( (fabs(trkQoverPt) > 1/minPtToReduceLayers_)  ?  (nLayers >= minStubLayers_)  :  (nLayers >= minStubLayers_ - 1) );
 } 

  //--- Filters returning filtered stubs based on input ones.
 
  // Produce a filtered collection of stubs from the input ones (on original track) that all have consistent rapidity
  std::vector<const Stub*> etaFilter ( const std::vector<const Stub*>& stubs, float trkQoverPt ) const;
  // Produce a filtered collection of stubs from the input ones (on original track) that all have consistent zR.
  std::vector<const Stub*> zTrkFilter (const std::vector<const Stub*>& stubs, float trkQoverPt );
  // Produce a filtered collection of stubs from the input ones (on original track)that are consistent with a straight line in r-z using tracklet algo.
  std::vector<const Stub*> seedFilter (const std::vector<const Stub*>& stubs, float trkQoverPt );

private:

  //=== Configuration parameters

  const Settings* settings_;

  float etaMinSector_; // rapidity range of this sector.
  float etaMaxSector_;
  float chosenRofZ_;   // Radius used to defined zTrkMinSector and zTrkMaxSector.
  float zTrkMinSector_; // corresponding range of this sector specified as z coordinate of track at given radius.
  float zTrkMaxSector_; 
  float phiCentreSector_; // phi coordinate of its centre.

  // Use filter in each HT cell using only stubs which have consistent rapidity?
  bool   useEtaFilter_;
  // Use filter in each HT cell using only stubs which have consistent zR
  bool   useZTrkFilter_;
  // Filter stubs in cell using a tracklet-like algorithm
  bool   useSeedFilter_;

  // Options for Ztrk filter
  float chosenRofZFilter_;

  // Options for Seed filter.
  float seedResolution_;
  bool  keepAllSeed_;

  // Number of seed combinations considered by the ZTrk Filter, for each input track.
  std::vector<unsigned int>  numZtrkSeedCombsPerTrk_;

  // Number of seed combinations considered by the Seed Filter, for each input track.
  std::vector<unsigned int>  numSeedCombsPerTrk_;
  std::vector<unsigned int>  numGoodSeedCombsPerTrk_;
  unsigned int maxSeedCombinations_;
  bool         zTrkSectorCheck_;

  // Min. numbers of layers with stubs required to produce track candidate.
  float  minStubLayers_; 
  // Reduce MinStubLayers by 1 for HT cells corresponding to track pt above this cut.
  // If this is set to >= 10000., this option is disabled.
  float  minPtToReduceLayers_;

  float beamWindowZ_; // Assumed length of beam spot in z.

  // Track (z0, tan_lambda) estimate from r-z filter if available.
  bool   estValid_;
  float  estZ0_;
  float  estTanLambda_;
};
#endif

