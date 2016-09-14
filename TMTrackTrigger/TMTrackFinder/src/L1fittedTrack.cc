#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"

using namespace std;

// Function for merging two tracks into a single track, used by by KillDupTracks.h for duplicate track removal.

L1fittedTrack L1fittedTrack::mergeTracks(const L1fittedTrack B) const {

  vector<const Stub*> aStubs=this->getStubs(), bStubs=B.getStubs();

  // Since a "set" contains no duplicates, add stubs to a set to eliminate the duplicates.
  set<const Stub*> mStubs;
  mStubs.insert(aStubs.begin(), aStubs.end());
  mStubs.insert(bStubs.begin(), bStubs.end());

  // Now copy the set back to the required vector for output.
  vector<const Stub*> mergedStubs;
  for (const Stub* s: mStubs) {
    mergedStubs.push_back(s);
  }

  // N.B. This defines the HT cell location as that of the first track, meaning that the merged tracks depends
  // on which track is first and which is second. This will make it hard to get identical results from hardware 
  // & software.
  return L1fittedTrack(settings_, l1track3D_, mergedStubs, 
		       qOverPt_, d0_, phi0_, z0_, tanLambda_, chi2_, nHelixParam_,
		       iPhiSec_, iEtaReg_, accepted_);
}
