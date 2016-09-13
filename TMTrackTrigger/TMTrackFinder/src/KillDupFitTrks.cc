#include "TMTrackTrigger/TMTrackFinder/interface/KillDupFitTrks.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <map>

//=== Make available cfg parameters & specify which algorithm is to be used for duplicate track removal.

void KillDupFitTrks::init(const Settings* settings, unsigned int dupTrkAlg)
{
  settings_ = settings;
  dupTrkAlg_ = dupTrkAlg;
  killDupTrks_.init(settings, dupTrkAlg); // Initialise duplicate removal algorithms that are common to all tracks.
}

//=== Eliminate duplicate tracks from the input collection, and so return a reduced list of tracks.

vector<L1fittedTrack> KillDupFitTrks::filter(const vector<L1fittedTrack>& vecTracks) const
{
  if (dupTrkAlg_ == 0) {

    // We are not running duplicate removal, so return original fitted track collection.
    // Note that this still includes tracks marked as "not accepted" by the fitter, but they are marked.
    // This allows the histogramming to plot also "not accepted" tracks.
    
    return vecTracks;

  } else {

    // We are running duplicate removal. It makes no sense to run this on tracks marked as "not accepted"
    // by the fitter, so remove them before proceeding.

    vector<L1fittedTrack> filtVecTracks;
    for (const L1fittedTrack& trk : vecTracks) {
      if (trk.accepted()) filtVecTracks.push_back(trk);
    }
	
    // Choose which algorithm to run, based on parameter dupTrkAlg_.
    switch (dupTrkAlg_) {
      // Run filters that only work on fitted tracks.
      case 50: return filterAlg50( filtVecTracks ); break;
      // Run filters that work on any type of track (l1track2d, l1track3d, l1fittedtrack). 
      default: return killDupTrks_.filter(filtVecTracks); 
    }
	
    // We should never end up here
    return filtVecTracks;
  }
}

//=== Duplicate removal algorithm designed to run after the track helix fit, which eliminates duplicates  
//=== simply by requiring that the fitted (q/Pt, phi0) of the track correspond to the same HT cell in 
//=== which the track was originally found by the HT.
//=== N.B. This code runs on tracks in a single sector. It could be extended to run on tracks in entire
//=== tracker by adding the track's sector number to memory "htCellUsed" below.


vector<L1fittedTrack> KillDupFitTrks::filterAlg50(const vector<L1fittedTrack>& tracks) const
{
  // Make a first pass through the tracks, doing initial identification of duplicate tracks.
  map<const L1fittedTrack*, bool> consistentMap; 
  set< pair<unsigned int, unsigned int> > htCellUsed;
  for (const L1fittedTrack& trk : tracks) {
    // Only consider tracks whose fitted helix parameters are in the same sector as the HT originally used to find the track.
    if (trk.consistentSector()) {
      // Check if this track's fitted (q/pt, phi0) helix parameters correspond to the same HT cell as the HT originally found the track in.
      bool consistentCell = trk.consistentHTcell();
      consistentMap[&trk] = consistentCell;  // Indicates if this track passes the HT cell consistency requirement.
      if (consistentCell) {
	// Memorize HT cell location corresponding to this track (identical for HT track & fitted track).
	htCellUsed.insert( trk.getL1track3D().getCellLocationRphi() );
      }
    }
  }

  // Making a second pass through the tracks, checking if any initially rejected should be rescued.
  vector<L1fittedTrack> tracksFiltered;
  for (const auto& mapPair : consistentMap) {
    const L1fittedTrack* trk         = mapPair.first;
    bool consistentCell  = mapPair.second;
    if (consistentCell) {
      // This track was selected by first pass, so keep it.
      tracksFiltered.push_back(*trk); 
    } else {
      // This track was rejected by first pass, so check if it should be rescued.
      // Get location in HT array corresponding to fitted track helix parameters.
      pair<unsigned int, unsigned int> htCell = trk->getCellLocationRphi();
      // If this HT cell was not already memorized, rescue this track, since it is probably not a duplicate,
      // but just a track whose fitted helix parameters are a bit wierd for some reason.
      if (std::count(htCellUsed.begin(), htCellUsed.end(), htCell) == 0) {
	tracksFiltered.push_back(*trk); // Rescue track.
      }
    }
  }

  return tracksFiltered;
}

