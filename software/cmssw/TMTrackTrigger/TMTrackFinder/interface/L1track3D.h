#ifndef __L1track3D_H__
#define __L1track3D_H__

#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1trackBase.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Utility.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"

#include <vector>
#include <utility>

using namespace std;

//=== L1 track candidate found in 3 dimensions.
//=== Gives access to all stubs on track and to its 3D helix parameters.
//=== Also calculates & gives access to associated truth particle (Tracking Particle) if any.

class L1track3D : public L1trackBase {

public:

  L1track3D(const Settings* settings, const vector<const Stub*>& stubs,
            pair<unsigned int, unsigned int> cellLocationRphi, pair<float, float> helixRphi,
            pair<unsigned int, unsigned int> cellLocationRz,   pair<float, float> helixRz) : 
    L1trackBase(),
    settings_(settings),
    stubs_(stubs), 
    cellLocationRphi_(cellLocationRphi), helixRphi_(helixRphi),
    cellLocationRz_  (cellLocationRz)  , helixRz_  (helixRz)
  {
    nLayers_   = Utility::countLayers(settings, stubs); // Count tracker layers these stubs are in
    matchedTP_ = Utility::matchingTP(settings, stubs, nMatchedLayers_, matchedStubs_); // Find associated truth particle & calculate info about match.
  }

  ~L1track3D() {}

  //--- Get information about the reconstructed track.

  // Get stubs on track candidate.
  const vector<const Stub*>&        getStubs()              const  {return stubs_;}  
  // Get number of stubs on track candidate.
  unsigned int                      getNumStubs()           const  {return stubs_.size();}
  // Get number of tracker layers these stubs are in.
  unsigned int                      getNumLayers()          const  {return nLayers_;}
  // Get cell location of track candidate in r-phi Hough Transform array in units of bin number.
  pair<unsigned int, unsigned int>  getCellLocationRphi()   const  {return cellLocationRphi_;}
  // The two conventionally agreed track helix parameters relevant in r-phi plane. i.e. (q/Pt, phi0)
  pair<float, float>                getHelixRphi()          const  {return helixRphi_;}
  // Get cell location of track candidate in r-z Hough Transform array in units of bin number.
  pair<unsigned int, unsigned int>  getCellLocationRz()     const  {return cellLocationRz_;}
  // The two conventionally agreed track helix parameters relevant in r-z plane. i.e. (z0, tan_lambda)
  pair<float, float>                getHelixRz()            const  {return helixRz_;}

  //--- User-friendlier access to the helix parameters. 

  float   qOverPt()    const  {return helixRphi_.first;}
  float   charge()     const  {return (this->qOverPt() > 0  ?  1  :  -1);} 
  float   invPt()      const  {return fabs(this->qOverPt());}
  float   pt()         const  {return 1./(1.0e-6 + this->invPt());} // includes protection against 1/pt = 0.
  float   d0()         const  {return 0.;} // Hough transform assumes d0 = 0.
  float   phi0()       const  {return helixRphi_.second;}
  float   z0()         const  {return helixRz_.first;}
  float   tanLambda()  const  {return helixRz_.second;}
  float   theta()      const  {return atan2(1., this->tanLambda());} // Use atan2 to ensure 0 < theta < pi.
  float   eta()        const  {return -log(tan(0.5*this->theta()));}

  // Phi and z coordinates at which track crosses "chosenR" values used by r-phi HT and rapidity sectors respectively.
  float   phiAtChosenR() const  {return reco::deltaPhi(this->phi0() - (settings_->invPtToDphi() * settings_->chosenRofPhi()) * this->qOverPt(),  0.);}
  float   zAtChosenR()   const  {return (this->z0() + (settings_->chosenRofZ()) * this->tanLambda());} // neglects transverse impact parameter & track curvature.

  // Comparitor for sorting tracks by q/Pt using std::sort().
  static bool qOverPtSortPredicate(const L1track3D& t1, const L1track3D t2) { return t1.getCellLocationRphi().first < t2.getCellLocationRphi().first; }

  //--- Get information about its association (if any) to a truth Tracking Particle.

  // Get best matching tracking particle (=nullptr if none).
  const TP*                  getMatchedTP()          const   {return matchedTP_;}
  // Get the matched stubs with this Tracking Particle
  const vector<const Stub*>& getMatchedStubs()       const   {return matchedStubs_;}
  // Get number of matched stubs with this Tracking Particle
  unsigned int               getNumMatchedStubs()    const   {return matchedStubs_.size();}
  // Get number of tracker layers with matched stubs with this Tracking Particle 
  unsigned int               getNumMatchedLayers()   const   {return nMatchedLayers_;}
  // Get purity of stubs on track candidate (i.e. fraction matching best Tracking Particle)
  float                      getPurity()             const   {return getNumMatchedStubs()/float(getNumStubs());}

  //--- Function for merging two tracks into a single track, used by by KillDupTracks.h for duplicate track removal.
  L1track3D mergeTracks(const L1track3D B) const;

private:

  //--- Configuration parameters
  const Settings*                    settings_; 

  //--- Information about the reconstructed track.
  vector<const Stub*>                stubs_;
  unsigned int                       nLayers_;
  pair<unsigned int, unsigned int>   cellLocationRphi_; 
  pair<float, float>                 helixRphi_; 
  pair<unsigned int, unsigned int>   cellLocationRz_; 
  pair<float, float>                 helixRz_; 

  //--- Information about its association (if any) to a truth Tracking Particle.
  const TP*                          matchedTP_;
  vector<const Stub*>                matchedStubs_;
  unsigned int                       nMatchedLayers_;
};
#endif
