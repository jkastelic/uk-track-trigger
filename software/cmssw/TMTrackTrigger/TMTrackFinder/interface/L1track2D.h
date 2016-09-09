#ifndef __L1track2D_H__
#define __L1track2D_H__

#include "FWCore/Utilities/interface/Exception.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1trackBase.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Utility.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"

#include <vector>
#include <utility>

using namespace std;

//=== L1 track cand found in 2 dimensions.
//=== Gives access to all stubs on track and to its 2D helix parameters.
//=== Also calculates & gives access to associated truth particle (Tracking Particle) if any.

class L1track2D : public L1trackBase {

public:

  // Give stubs on track, its cell location inside HT arraym its 2D helix parameters, and 
  // indicate if it was found by an r-phi or r-z HT. 
  L1track2D(const Settings* settings, const vector<const Stub*>& stubs, 
	    pair<unsigned int, unsigned int> cellLocation, pair<float, float> helix2D, bool isRphi) : 
            L1trackBase(),
	    settings_(settings),
	    stubs_(stubs), 
            cellLocation_(cellLocation),
            helix2D_(helix2D),
	    isRphi_(isRphi),
	    estValid_(false),
	    estZ0_(0.),
	    estTanLambda_(0.)
  {
    nLayers_   = Utility::countLayers(settings, stubs); // Count tracker layers these stubs are in
    matchedTP_ = Utility::matchingTP(settings, stubs, nMatchedLayers_, matchedStubs_); // Find associated truth particle & calculate info about match.
  }

  ~L1track2D() {}

  //--- Get information about the reconstructed track.

  // Get stubs on track candidate.
  const vector<const Stub*>&        getStubs()              const  {return stubs_;}  
  // Get number of stubs on track candidate.
  unsigned int                      getNumStubs()           const  {return stubs_.size();}
  // Get number of tracker layers these stubs are in.
  unsigned int                      getNumLayers()          const  {return nLayers_;}
  // Get cell location of track candidate in Hough Transform array in units of bin number.
  pair<unsigned int, unsigned int>  getCellLocation()       const  {return cellLocation_;}
  // The two conventionally agreed track helix parameters relevant in this 2D plane.
  // i.e. (q/Pt, phi0) or (z0, tan_lambda).
  pair<float, float>                getHelix2D()            const  {return helix2D_;}
  // Was it found by an r-phi or r-z Hough transform?
  bool                              isRphiTrk()             const  {return isRphi_;}

  //--- User-friendlier access to the helix parameters obtained from track location inside HT array.

  float   qOverPt()    const  {if (! isRphi_) throw cms::Exception("L1track2D ERROR: You asked for r-phi helix params of r-z track."); return helix2D_.first;}
  float   phi0()       const  {if (! isRphi_) throw cms::Exception("L1track2D ERROR: You asked for r-phi helix params of r-z track."); return helix2D_.second;}
  float   z0()         const  {if (isRphi_)   throw cms::Exception("L1track2D ERROR: You asked for r-z helix params of r-phi track."); return helix2D_.first;}
  float   tanLambda()  const  {if (isRphi_)   throw cms::Exception("L1track2D ERROR: You asked for r-z helix params of r-phi track."); return helix2D_.second;}

  // Comparitor for sorting tracks by q/Pt using std::sort().
  static bool qOverPtSortPredicate(const L1track2D& t1, const L1track2D t2) { return t1.getCellLocationRphi().first < t2.getCellLocationRphi().first; }

  //--- In the case of tracks found by the r-phi HT, a rough estimate of the (z0, tan_lambda) may be provided by any r-z
  //--- track filter run after the r-phi HT. These two functions give set/get access to these.
  //--- The "get" function returns a boolean indicating if an estimate exists (i.e. "set" has been called).

  void setTrkEstZ0andTanLam(float  estZ0, float  estTanLambda) {
    if (! isRphi_) throw cms::Exception("L1track2D ERROR: You provided estimated r-z helix paramers for an r-z track.");
    estZ0_ = estZ0; estTanLambda_ = estTanLambda; estValid_ = true;
  } 
  bool getTrkEstZ0andTanLam(float& estZ0, float& estTanLambda) const {
    if (! isRphi_) throw cms::Exception("L1track2D ERROR: You tried to retrieve estimated r-z helix paramers for an r-z track.");
    estZ0 = estZ0_; estTanLambda = estTanLambda_; return estValid_;
  }

  //--- User-friendlier access to the cell locations of the track candidate in the r-phi and r-z Hough transform arrays in units of bin number.

  pair<unsigned int, unsigned int>  getCellLocationRphi() const {if (! isRphi_) throw cms::Exception("L1track2D ERROR: You asked for r-phi HT cell of r-z track."); return cellLocation_;}
  pair<unsigned int, unsigned int>  getCellLocationRz()   const {if (  isRphi_) throw cms::Exception("L1track2D ERROR: You asked for r-z HT cell of r-phi track."); return cellLocation_;}

  //--- Get information about its association (if any) to a truth Tracking Particle.

  // Get matching tracking particle (=nullptr if none).
  const TP*                  getMatchedTP()          const   {return matchedTP_;}
  // Get the matched stubs.
  const vector<const Stub*>& getMatchedStubs()       const   {return matchedStubs_;}
  // Get number of matched stubs.
  unsigned int               getNumMatchedStubs()    const   {return matchedStubs_.size();}
  // Get number of tracker layers with matched stubs.
  unsigned int               getNumMatchedLayers()   const   {return nMatchedLayers_;}

  //--- Function for merging two tracks into a single track, used by by KillDupTracks.h for duplicate track removal.
  L1track2D mergeTracks(const L1track2D B) const;

private:

  //--- Configuration parameters
  const Settings*                    settings_; 

  //--- Information about the reconstructed track from Hough transform.
  vector<const Stub*>                stubs_;
  unsigned int                       nLayers_;
  pair<unsigned int, unsigned int>   cellLocation_; 
  pair<float, float>                 helix2D_; 
  bool                               isRphi_;

  //--- Rough estimate of r-z track parameters from r-z filter, which may be present in case of r-phi Hough transform
  bool  estValid_;
  float estZ0_;
  float estTanLambda_;

  //--- Information about its association (if any) to a truth Tracking Particle.  
  const TP*                          matchedTP_;
  vector<const Stub*>                matchedStubs_;
  unsigned int                       nMatchedLayers_;
};
#endif
