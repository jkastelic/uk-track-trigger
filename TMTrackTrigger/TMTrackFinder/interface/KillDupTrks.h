#ifndef __KILLDUPTRKS_H__
#define __KILLDUPTRKS_H__

#include <cstddef>
#include <vector>
#include <algorithm>
#include <functional>
#include <utility>
#include <iostream>
#include <type_traits>

#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <gsl/gsl_fit.h>

using namespace std;

class L1trackBase;
class L1track2D;
class L1track3D;
class L1fittedTrack;

/**
*  Kill duplicate reconstructed tracks.
*  e.g. Those sharing many hits in common.
*  
*  Currently this is intended to run only on tracks found within a single (eta,phi) sector.
*  It includes a naive algorithms from Ian (dupTrkAlg = 1) & more sophisticated ones from Ivan (dupTrkAlg > 1).
*  The class is implemented inside TMTrackTrigger/TMTrackFinder/interface/KillDupTrks.icc
*  
*  The template class "T" can be any class inheriting from L1trackBase.
* 
*  -------------------------------------------------------------------------------------------
*   GENERAL INFO ABOUT THE FILTER ALGORITHMS DEFINED IN THE CLASS.
*   Some of these algorithms are designed to work on r-phi L1track2D tracks, and some on r-z 
*   L1track2D tracks. Others work on L1tracks3D.
*  -------------------------------------------------------------------------------------------
*/
template <class T> class KillDupTrks {

public:

  KillDupTrks()
	{
    // Check that classed used as template "T" inherits from class L1trackBase.
    static_assert(std::is_base_of<L1trackBase, T>::value, "KillDupTrks ERROR: You instantiated this with a template class not inheriting from L1trackBase!");
  }

  ~KillDupTrks() {}

  /**
  *  Make available cfg parameters & specify which algorithm is to be used for duplicate track removal.
  */
  void init(const Settings* settings, unsigned int dupTrkAlg);

  /**
  *  Eliminate duplicate tracks from the input collection, and so return a reduced list of tracks.
  */
  vector<T> filter(const vector<T>& vecTracks) const;

private:


  /**
  *  A specific algorithm for filtering duplicate tracks.
  *  Selects only a single track, based on the number of stubs & number of layers with stubs on the track.
  */
  vector<T> filterAlg1(const vector<T>& vecTracks) const;

  /**
  *  An algorithm to remove candidates with exactly the same stubs as another
  *  Based on Stub index() -- assumes they are ordered!  idr 9/7/15
  */
  vector<T> filterAlg2(const vector<T>& vecTracks) const;

  /**
  *  A specific algorithm for filtering duplicates
  *  Pairwise candidate comparison, removes tracks with fewer than N independent stubs
  *  Implementing OSU algorithm, keep tracks with N or more unique stubs (default 3)
  */
  vector<T> filterAlg3(const vector<T>& vecTracks) const;

  /**
  *  A specific algorithm for filtering duplicates
  *  Cut on ChiSq of linear fit in RZ -- experimental, didn't work well
  *  Filter on reduced ChiSq of a linear fit in RZ
  */
  vector<T> filterAlg4(const vector<T>& vecTracks) const;

  /**
  *  A specific algorithm for filtering duplicates
  *  Pairwise candidate comparison, if two have at least N common stubs in N layers, keep (smaller) larger
  *  Implementing "inverse" OSU algorithm, check for stubs in common,
  *  keep {smallest} longest! candidates if common stubs in N or more layers (default 5 at present)
  */
  vector<T> filterAlg5(const vector<T>& vecTracks) const;

  /**
  *  A specific algorithm for filtering duplicates
  *  Pairwise candidate comparison, if two have at least N common stubs in N layers, keep one with smaller RZ/ZR red chisq
  *  Implementing "inverse" OSU algorithm, check for stubs in common,
  *  keep smallest candidates if common stubs in N or more layers (default 5 at present)
  *  Try keeping track with best RZ/ZR reduced chi-square
  */
  vector<T> filterAlg6(const vector<T>& vecTracks) const;

  /**
  *  A specific algorithm for filtering duplicates
  *  Pairwise candidate comparison, if two have at least N common stubs in N layers, keep one with best "quality"
  *  Implementing "inverse" OSU algorithm, check for stubs in common,
  *  keep "best"((c) IanT 2015) candidates if common stubs in N or more layers (default 5 at present)
  */
  vector<T> filterAlg7(const vector<T>& vecTracks) const;

  /**
  *  Implementing "inverse" OSU algorithm, check for stubs in common,
  *  keep largest candidates if common stubs in N or more layers (default 5 at present), both if equal
  *  Implementing "inverse" OSU algorithm, check for stubs in common,
  *  keep largest candidates if common stubs in N or more layers (default 5 at present), both if equal
  */
  vector<T> filterAlg8(const vector<T>& vecTracks) const;

  /**
  *  Implementing "inverse" OSU algorithm, check for stubs in common,
  *  keep smallest candidates if common stubs in N or more layers (default 5 at present),
  *  Didn't work, back to Alg8 for present...
  *  later add keep least stubs if equal
  *  Implementing "inverse" OSU algorithm, check for stubs in common,
  *  Try keeping _smallest_; Nope, didn't work, back to original Alg8
  *  (later -- if equal, keep least no of stubs), if still equal discard one (otherwise dupes not removed)
  */
  vector<T> filterAlg9(const vector<T>& vecTracks) const;

  /**
  *  Try just removing r-phi candidates in adjacent cells as dupes are mostly adjacent
  */
  vector<T> filterAlg10(const vector<T>& vecTracks) const;

  /**
  *  Try just removing r-phi candidates in adjacent cells *with same number of stubs* as dupes are mostly adjacent
  */
  vector<T> filterAlg11(const vector<T>& vecTracks) const;

  /**
  *  Try just removing r-phi candidates in adjacent cells in X (with same number of stubs) as dupes are mostly adjacent
  */
  vector<T> filterAlg12(const vector<T>& vecTracks) const;

  /**
  *  Try just merging r-phi candidates in adjacent cells in X
  */
  vector<T> filterAlg13(const vector<T>& vecTracks) const;

  /**
  *  Try merging r-phi candidates in all adjacent cells
  */
  vector<T> filterAlg14(const vector<T>& vecTracks) const;

  /**
  *  Implementing "inverse" OSU algorithm, check for stubs in common,
  *  keep largest candidates if common stubs in N or more layers (default 5 at present)
  *  limit range of qOverPt scanned for duplicates
  */
  vector<T> filterAlg15(const vector<T>& vecTracks) const;

  /**
  *  Implementing "inverse" OSU algorithm, check for stubs in common,
  *  merge candidates if common stubs in N or more layers (default 5 at present)
  */
  vector<T> filterAlg16(const vector<T>& vecTracks) const;

  /**
   * Duplicate removal algorithm which merges candidates within cuts on the differences
   * between the four helix parameters
   */
  vector<T> filterAlg17(const vector<T>& vecTracks) const;

  /**
   * Duplicate removal algorithm which merges candidates in adjacent cells within cuts on the differences
   * between the dz0 and dtanLambda helix parameters
   */
  vector<T> filterAlg18(const vector<T>& vecTracks) const;

 /**
   * Duplicate removal algorithm which deletees candidates in adjacent cells within cuts on the differences
   * between the dz0 and dtanLambda helix parameters if worse "quality" (trying number of stubs...)(now layers...)
   */
  vector<T> filterAlg19(const vector<T>& vecTracks) const;

  /**
  *  Use a hash to identify whether stubs are identical between each other. The hash maps 
	*  to a space addressable with 16 bits (i.e. STL C++ hash is used modulo 2**16).
  */
  vector<T> filterAlg100(const vector<T>& vecTracks) const;


  /**
  *  Prints out a consistently formatted formatted report of killed duplicate track
  */
  void printKill(unsigned alg, unsigned dup, unsigned cand, T dupTrack, T candTrack) const;

  /**
  * Tests if cells are adjacent in q/pT
  */
  bool isNextQoverPt(std::pair<unsigned int, unsigned int> i, std::pair<unsigned int, unsigned int> j) const;

  /**
  * Tests if cells are adjacent in q/pT AND phi0
  */
  bool isAdjacentCell(std::pair<unsigned int, unsigned int> i, std::pair<unsigned int, unsigned int> j) const;

private:

  const Settings *settings_; // Configuration parameters.

  bool  enableMerge2x2_;
  float minInvPtToMerge2x2_;
  float maxAbsQoverPtAxis_;       // Max. |q/Pt| covered by  HT array.
  unsigned int nBinsQoverPtAxis_; // Number of bins in HT array in q/Pt.
  unsigned int minFineBin_;       // Start of finer binning
  unsigned int maxFineBin_;       // End of finer binning

  unsigned int dupTrkAlg_; // Specifies choice of algorithm for duplicate track removal.
  unsigned int dupTrkMinIndependent_; // Minimum number of independent stubs to retain track in OSU alg
  unsigned int dupTrkMinCommonHitsLayers_;  // Min no of matched stubs & layers to keep smaller cand
  double dupTrkChiSqCut_;  // Value for reduced ChiSq cut in Alg4
  float dupMaxQOverPtScan_;  //Cutoff for qOverPt difference in Algorithm 15
  float dupMaxPhi0Scan_;  //Cutoff for phi0 in Alg 15 etc.
  float dupMaxZ0Scan_;  //Cutoff for z0 in Alg 15 etc.
  float dupMaxTanLambdaScan_; //Cutoff for tan lambda in Alg 15 etc.
};

//=== Include file which implements all the functions in the above class.
#include "TMTrackTrigger/TMTrackFinder/interface/KillDupTrks.icc"

#endif

