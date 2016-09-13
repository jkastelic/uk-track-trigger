#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/DigitalStub.h"

#include "DataFormats/Math/interface/deltaPhi.h"

//--- Note that in all that follows, anything related to dphi & rho is irrelevant when using the "daisy chain"
//--- firmware. The digitization is done below for these variables, but is never used for daisu chain.

//=== Note configuration parameters.

DigitalStub::DigitalStub(const Settings* settings) :

  // To checkn that DigitalStub is correctly initialized.
  ranInit_(false),
  ranMakeGPinput_(false),
  ranMakeHTinput_(false),

  // Digitization configuration parameters
  iFirmwareType_ (settings->firmwareType()),  // Firmware type (systolic, Thomas ...)
  phiSectorBits_ (settings->phiSectorBits()), // No. of bits to store phi sector number
  //--- Parameters available in HT board.
  phiSBits_      (settings->phiSBits()),      // No. of bits to store phiS coord.
  phiSRange_     (settings->phiSRange()),     // Range of phiS coord. in radians.
  rtBits_        (settings->rtBits()),        // No. of bits to store rT coord.
  rtRange_       (settings->rtRange()),       // Range of rT coord. in cm.
  zBits_         (settings->zBits()),         // No. of bits to store z coord.
  zRange_        (settings->zRange()),        // Range of z coord in cm.
  dPhiBits_      (settings->dPhiBits()),      // No. of bits to store Delta(phi).
  dPhiRange_     (settings->dPhiRange()),     // Range of Delta(phi) in radians.
  rhoBits_       (settings->rhoBits()),       // No. of bits to store rho parameter.
  rhoRange_      (settings->rhoRange()),      // Range of rho parameter.
  //--- Parameters available in GP board (excluding any in common with HT specified above).
  phiOBits_      (settings->phiOBits()),      // No. of bits to store phiO parameter.
  phiORange_     (settings->phiORange()),     // Range of phiO parameter
  bendBits_      (settings->bendBits()),      // No. of bits to store stub bend.

  // Note if using reduced layer ID, so tracker layer can be encoded in 3 bits.
  reduceLayerID_ (settings->reduceLayerID()),

  // Number of phi sectors and phi octants.
  numPhiSectors_ (settings->numPhiSectors()),
  numPhiOctants_ (8),
  // Phi sector and phi octant width (radians)
  phiSectorWidth_(2.*M_PI / float(numPhiSectors_)), 
  phiOctantWidth_(2.*M_PI / float(numPhiOctants_)), 
  // Radius from beamline with respect to which stub r coord. is measured.
  chosenRofPhi_  (settings->chosenRofPhi()),

  // Number of q/Pt bins in Hough  transform array.
  nbinsPt_       ((int) settings->houghNbinsPt())
{
  // Calculate multipliers to digitize the floating point numbers.
  phiSMult_ = pow(2, phiSBits_)/phiSRange_;
  rtMult_   = pow(2, rtBits_  )/rtRange_;
  zMult_    = pow(2, zBits_   )/zRange_;
  dPhiMult_ = pow(2, dPhiBits_)/dPhiRange_;
  rhoMult_  = pow(2, rhoBits_ )/rhoRange_;
  phiOMult_ = pow(2, phiOBits_)/phiORange_;

  bendMult_ = 4.; // No precision lost by digitization, since original bend (after encoding) has steps of 0.25*pitch
  bendRange_ = pow(2, bendBits_)/bendMult_; 
}

//=== Initialize stub with original, floating point stub coords, stub bend angle, and rho parameter,
//=== range of m bin (= q/Pt bin) values allowed by bend filter, 
//=== normal & "reduced" tracker layer of stub, stub bend, and pitch & seperation of module.

void DigitalStub::init(float phi_orig, float r_orig, float z_orig, float dphi_orig, float rho_orig,
		       unsigned int min_qOverPt_bin_orig, unsigned int max_qOverPt_bin_orig, 
		       unsigned int layerID, unsigned int layerIDreduced, float bend_orig,
		       float pitch, float sep) {

  ranInit_ = true; // Note we ran init().
  // Variables in HT.
  phi_orig_             = phi_orig; 
  r_orig_               = r_orig;
  z_orig_               = z_orig;
  dphi_orig_            = dphi_orig;
  rho_orig_             = rho_orig;
  min_qOverPt_bin_orig_ = min_qOverPt_bin_orig;
  max_qOverPt_bin_orig_ = max_qOverPt_bin_orig;
  layerID_              = layerID;
  layerIDreduced_       = layerIDreduced;
  // Variables exclusively in GP.
  bend_orig_            = bend_orig;

  // Calculate unique module type ID, allowing pitch/sep of module to be determined.

  bool barrel = (layerID < 10);
  moduleType_ = 999;
  const vector<float> pitchVsType  = {0.01, 0.01, 0.01, 0.009, 0.009, 0.009};
  const vector<float> sepVsType    = {0.26, 0.16, 0.4, 0.18, 0.18, 0.4};
  const vector<bool>  barrelVsType = {true, true, false, true, false, false};
  if (pitchVsType.size() != sepVsType.size()) throw cms::Exception("DigitalStub: module type array size wrong");
  const float tol = 0.001; // Tolerance
  for (unsigned int i = 0; i < pitchVsType.size(); i++) {
    if (fabs(pitch - pitchVsType[i]) < tol && fabs(sep - sepVsType[i]) < tol && barrel == barrelVsType[i]) {
      moduleType_ = i;
    }
  }
  if (moduleType_ == 999) throw cms::Exception("DigitalStub: unknown module type")<<"pitch="<<pitch<<" separation="<<sep<<" barrel="<<barrel<<endl;
}

//=== Digitize stub for input to Geographic Processor, with stub phi coord. measured relative to phi octant that contains specified phi sector.

void DigitalStub::makeGPinput(unsigned int iPhiSec) {

  if (! ranInit_) throw cms::Exception("DigitalStub:: You forgot to call init() before make()!");

  unsigned int iPhiOct = floor(iPhiSec*numPhiOctants_/numPhiSectors_); // Find octant corresponding to this sector.

  // If this stub was already digitized, we don't have to redo all the work again. Save CPU.
  if (ranMakeGPinput_) {
    if (iPhiOct == iDigi_Octant_) {
      return; // Work already done.
    } else {
      this->quickMakeGPinput(iPhiSec);
      return;
    }
  }

  ranMakeGPinput_ = true; // Note we ran make().

  //--- Shift axes of coords. if required.

  // r coordinate relative to specified point.
  rt_orig_ = r_orig_ - chosenRofPhi_;

  // Phi coord. of stub relative to centre of octant.
  double phiOctantCentre = phiOctantWidth_ * (0.5 + double(iPhiOct)) - M_PI;
  phiO_orig_ = reco::deltaPhi(phi_orig_, phiOctantCentre);

  //--- Digitize variables used exclusively in GP.
  iDigi_Octant_ = iPhiOct;
  iDigi_PhiO_   = floor(phiO_orig_*phiOMult_);
  iDigi_Bend_   = floor(bend_orig_*bendMult_); 
  //--- Digitize variables used in both GP & HT.
  iDigi_Rt_     = floor(rt_orig_*rtMult_);
  iDigi_Z_      = floor(z_orig_*zMult_);
  
  //--- Determine floating point stub coords. from digitized numbers (so with degraded resolution).
  //--- First for variables used exclusively in GP.
  phiO_      = (iDigi_PhiO_ + 0.5)/phiOMult_;
  bend_      = (iDigi_Bend_      )/bendMult_; // Different eqn. as bend is known to be half-integer.
  phi_       = reco::deltaPhi(phiO_, -phiOctantCentre);
  //--- Then for variables used in both GP & HT.
  rt_        = (iDigi_Rt_   + 0.5)/rtMult_;  
  r_         = rt_ + chosenRofPhi_;
  z_         = (iDigi_Z_    + 0.5)/zMult_; 
}

//=== Digitize stub for input to Hough transform, with stub phi coord. measured relative to specified phi sector.

void DigitalStub::makeHTinput(unsigned int iPhiSec) {

  if (! ranInit_) throw cms::Exception("DigitalStub:: You forgot to call init() before make()!");

  // Digitize for GP input if not already done, since some variables are shared by GP & HT.
  this->makeGPinput(iPhiSec);

  // If this stub was already digitized, we don't have to redo all the work again. Save CPU.
  if (ranMakeHTinput_) {
    if (iPhiSec == iDigi_PhiSec_) {
      return; // Work already done.
    } else {
      this->quickMakeHTinput(iPhiSec);
      return;
    }
  }

  ranMakeHTinput_ = true; // Note we ran makeHTinput().

  //--- Shift axes of coords. if required.

  // Centre of this sector in phi
  double phiSectorCentre = phiSectorWidth_ * (0.5 + double(iPhiSec)) - M_PI; 
  // Point in sector from which stub phiS should be measured.
  double phiSectorRef = phiSectorCentre;
  // Systolic array firmware measures might measure phiS from start of sector, not centre.
  if (iFirmwareType_ == 9) phiSectorRef = phiSectorCentre - phiSectorWidth_*0.5;

  // Phi coord. of stub relative to centre of sector.
  phiS_orig_ = reco::deltaPhi(phi_orig_, phiSectorRef); 

  //--- Digitize variables used exclusively by HT.
  iDigi_PhiSec_ = iPhiSec;
  iDigi_PhiS_   = floor(phiS_orig_*phiSMult_);
  iDigi_Dphi_   = floor(dphi_orig_*dPhiMult_);
  iDigi_Rho_    = floor(rho_orig_*rhoMult_);
  // Don't bother digitising here variables used by both GP & HT, as makeGPinput() will already have digitized them.
  
  //--- Determine floating point stub coords. from digitized numbers (so with degraded resolution).
  //--- First for variables used exclusively by HT
  phiS_      = (iDigi_PhiS_ + 0.5)/phiSMult_;
  phi_       = reco::deltaPhi(phiS_, -phiSectorRef); // N.B. phi_ measured w.r.t sector here, but w.r.t. octant in makeGPinput()
  dphi_      = (iDigi_Dphi_ + 0.5)/dPhiMult_;
  rho_       = (iDigi_Rho_  + 0.5)/rhoMult_;
  // Don't bother with  variables used by both GP & HT, as makeGPinput() will already have already calculated them.

  //--- Do next two checks here rather than in makeGPinput(), as the latter may be called for stubs in the wrong octant, so out of range.

  // Check that stub coords. are within assumed digitization range.
  this->checkInRange();

  // Check that digitization followed by undigitization doesn't change results too much.
  this->checkAccuracy();

  // Adjust m-bin range, as hardware counts q/Pt bins in HT array using a signed integer in a symmetric range about zero.
  const int min_array_bin = (nbinsPt_%2 == 0)  ?  -(nbinsPt_/2)      :  -(nbinsPt_ - 1)/2;
  m_min_ = min_qOverPt_bin_orig_ + min_array_bin;
  m_max_ = max_qOverPt_bin_orig_ + min_array_bin;

  //--- Produce tracker layer identifier, encoded as it is sent along the optical link.
  
  if (reduceLayerID_) {
    // Firmware is using "reduced" layer ID, which can be packed into 3 bits in range 1-7.
    iDigi_LayerID_ = layerIDreduced_;

  } else {
    // Firmware is using normal layer ID, which needs more than 3 bits to store it.
    // Encode barrel layers as 0 to 5.
    iDigi_LayerID_ = layerID_ - 1;
    // Endcode endcap layers as 6 to 10, not bothering to distinguish the two endcaps.
    if        (iDigi_LayerID_ == 10 || iDigi_LayerID_ == 20) {
      iDigi_LayerID_ = 6;
    } else if (iDigi_LayerID_ == 11 || iDigi_LayerID_ == 21) {
      iDigi_LayerID_ = 7;
    } else if (iDigi_LayerID_ == 12 || iDigi_LayerID_ == 22) {
      iDigi_LayerID_ = 8;
    } else if (iDigi_LayerID_ == 13 || iDigi_LayerID_ == 23) {
      iDigi_LayerID_ = 9;
    } else if (iDigi_LayerID_ == 14 || iDigi_LayerID_ == 24) {
      iDigi_LayerID_ = 10;
    }
  }
}

//=== Redigitize stub for input to Geographic Processor, if it was previously digitized for a different phi sector.

void DigitalStub::quickMakeGPinput(int iPhiSec) {

  //--- Shift axes of coords. if required.

  // Phi coord. of stub relative to centre of octant.
  unsigned int iPhiOct = floor(iPhiSec*numPhiOctants_/numPhiSectors_);
  double phiOctantCentre = phiOctantWidth_ * (0.5 + double(iPhiOct)) - M_PI;
  phiO_orig_ = reco::deltaPhi(phi_orig_, phiOctantCentre);

  //--- Digitize variables used exclusively in GP.
  iDigi_Octant_ = iPhiOct;
  iDigi_PhiO_   = floor(phiO_orig_*phiOMult_);
  
  //--- Determine floating point stub coords. from digitized numbers (so with degraded resolution).
  //--- Variables used exclusively in GP.
  phiO_      = (iDigi_PhiO_ + 0.5)/phiOMult_;
  phi_       = reco::deltaPhi(phiO_, -phiOctantCentre);
}

//=== Redigitize stub for input to Hough transform, if it was previously digitized for a different phi sector.

void DigitalStub::quickMakeHTinput(int iPhiSec) {

  //--- Shift axes of coords. if required.

  // Centre of this sector in phi
  double phiSectorCentre = phiSectorWidth_ * (0.5 + double(iPhiSec)) - M_PI; 
  // Point in sector from which stub phiS should be measured.
  double phiSectorRef = phiSectorCentre;
  // Systolic array firmware measures might measure phiS from start of sector, not centre.
  if (iFirmwareType_ == 9) phiSectorRef = phiSectorCentre - phiSectorWidth_*0.5;

  // Phi coord. of stub relative to centre of sector.
  phiS_orig_ = reco::deltaPhi(phi_orig_, phiSectorRef); 

  // Check that stub coords. are within assumed digitization range.
  this->checkInRange();

  //--- Digitize variables used in HT.
  iDigi_PhiSec_ = iPhiSec;
  iDigi_PhiS_   = floor(phiS_orig_*phiSMult_);
  
  //--- Determine floating point stub coords. from digitized numbers (so with degraded resolution).
  //--- First for variables used in HT.
  phiS_      = (iDigi_PhiS_ + 0.5)/phiSMult_;
  phi_       = reco::deltaPhi(phiS_, -phiSectorRef);
}

//=== Check that stub coords. are within assumed digitization range.

void DigitalStub::checkInRange() const {
  // All ranges are centred at zero, except for rho, which is +ve-definate.
  if (fabs(phiS_orig_) >= 0.5*phiSRange_)   throw cms::Exception("DigitalStub: Stub phiS is out of assumed digitization range.")<<" |phiS| = " <<fabs(phiS_orig_) <<" > "<<0.5*phiSRange_<<endl;  
  if (fabs(rt_orig_)   >= 0.5*rtRange_)     throw cms::Exception("DigitalStub: Stub rT is out of assumed digitization range.")  <<" |rt| = "   <<fabs(rt_orig_)   <<" > "<<0.5*rtRange_  <<endl;  
  if (fabs(z_orig_)    >= 0.5*zRange_)      throw cms::Exception("DigitalStub: Stub z is out of assumed digitization range.")   <<" |z| = "    <<fabs(z_orig_)   <<" > "<<0.5*zRange_  <<endl;  
  if (fabs(dphi_orig_)    >= 0.5*dPhiRange_) throw cms::Exception("DigitalStub: Stub dphi is out of assumed digitization range.")<<" |dphi| = " <<fabs(dphi_orig_)<<" > "<<0.5*dPhiRange_  <<endl;  
  if (rho_orig_ <= 0. || rho_orig_ >= rhoRange_)   throw cms::Exception("DigitalStub: Stub rho is out of assumed digitization range.") <<" rho = "    <<rho_orig_ <<" is < 0 or > "<<rhoRange_<<endl;  
  if (fabs(phiO_orig_)   >= 0.5*phiORange_) throw cms::Exception("DigitalStub: Stub phiO is out of assumed digitization range.")<<" |dphi| = " <<fabs(phiO_orig_)<<" > "<<0.5*phiORange_  <<endl;  
  if (fabs(bend_orig_)   >= 0.5*bendRange_) throw cms::Exception("DigitalStub: Stub bend is out of assumed digitization range.")<<" |bend| = " <<fabs(bend_orig_)<<" > "<<0.5*bendRange_  <<endl;  
}

//=== Check that digitisation followed by undigitisation doesn't change significantly the stub coordinates.

void DigitalStub::checkAccuracy() const {
  float TA = reco::deltaPhi(phi_, phi_orig_);
  float TB = r_    - r_orig_;
  float TC = z_    - z_orig_;
  float TD = dphi_ - dphi_orig_;
  float TE = rho_  - rho_orig_;
  float TF = phiO_ - phiO_orig_;
  float TG = bend_ - bend_orig_;

  if (fabs(TA) > 0.001 || fabs(TB) > 0.3 || fabs(TC) > 0.2 || fabs(TD) > 0.005 || fabs(TE) > 0.005 || fabs(TF) > 0.005 || fabs(TG) > 0.01) {
    cout<<"STUB DIGI MESS UP "<<TA<<" "<<TB<<" "<<TC<<" "<<TD<<" "<<TE<<" "<<TF<<" "<<TG<<endl;
  }
}
