#include "TMTrackTrigger/TMTrackFinder/interface/L1track2D.h"

// Function for merging two tracks into a single track, used by by KillDupTracks.h for duplicate track removal.

L1track2D L1track2D::mergeTracks(const L1track2D B) const
{
	using namespace std;
	
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
  return L1track2D(settings_, mergedStubs, this->getCellLocation(), this->getHelix2D(), this->isRphiTrk());
}
