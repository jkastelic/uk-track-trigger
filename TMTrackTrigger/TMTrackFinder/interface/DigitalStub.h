#ifndef __DIGITALSTUB_H__
#define __DIGITALSTUB_H__

#include <cmath>

#include "FWCore/Utilities/interface/Exception.h"

using namespace std;

class Settings;

//=== Used to digitize stubs both for input to GP and for input to HT.
//=== N.B. After constructing an object of type DigitalStub, you must call functions 
//=== init() and make() before you try using any of the other functions to access digitized stub info.

//=== WARNING: Not all variables available in the GP are available inside the HT or visa-versa,
//=== so think about the hardware when calling the functions below.

class DigitalStub {

public:

  // Note configuration parameters.
  DigitalStub(const Settings* settings);

  ~DigitalStub(){}

  // Initialize stub with original, floating point stub coords, stub bend angle, and rho parameter,
  // range of m bin (= q/Pt bin) values allowed by bend filter, 
  // normal & "reduced" tracker layer of stub, stub bend, and pitch & seperation of module.
  void init(float phi_orig, float r_orig, float z_orig, float dphi_orig, float rho_orig, 
	    unsigned int min_qOverPt_bin_orig, unsigned int max_qOverPt_bin_orig, 
	    unsigned int layerID, unsigned int layerIDreduced, float bend_orig,
	    float pitch, float sep);

  // Digitize stub for input to Geographic Processor, with stub phi coord. measured relative to phi octant that contains specified phi sector.
  void makeGPinput(unsigned int iPhiSec);

  // Digitize stub for input to Hough transform, with stub phi coord. measured relative to specified phi sector.
  // (Also still allows digital data for input to GP to be accessed).
  void makeHTinput(unsigned int iPhiSec);

  //--- Note that the variables related to "dphi" and "rho" are not digitized when using the "daisy chain"
  //--- firmware, so no attempt should be made to access them in this case. They are replaced by the
  //--- m_min and m_max variables.

  //--- The functions below return variables post-digitization.
  //--- Do not call any of the functions below, unless you have already called init() and make()!

  // Digits corresponding to stub coords.
  // %%% Those common to GP & HT input.
  int          iDigi_Rt()      const {this->okGP(); return iDigi_Rt_;}     // r coord. relative to chosen radius
  int          iDigi_Z()       const {this->okGP(); return iDigi_Z_;}      // z coord.
  // %%% Those exclusively input to HT.
  unsigned int iDigi_PhiSec()  const {this->okHT(); return iDigi_PhiSec_;} // phi sector number
  int          iDigi_PhiS()    const {this->okHT(); return iDigi_PhiS_;}   // phi coord. relative to sector
  int          iDigi_Dphi()    const {this->okHT(); this->valid(); return iDigi_Dphi_;}   // Delta(phi) = bend angle of stub
  unsigned int iDigi_Rho()     const {this->okHT(); this->valid(); return iDigi_Rho_;}    // Rho
  // %%% Those exclusively input to GP.
  unsigned int moduleType()    const {this->okGP(); return moduleType_;}   // module type ID (gives pitch/spacing)
  unsigned int iDigi_Octant()  const {this->okGP(); return iDigi_Octant_;} // phi octant number
  int          iDigi_PhiO()    const {this->okGP(); return iDigi_PhiO_;}   // phi coord. relative to octant
  int          iDigi_Bend()    const {this->okGP(); return iDigi_Bend_;}   // stub bend

  // Floating point stub coords derived from digitized info (so with degraded resolution).
  // %%% Those common to GP & HT input.
  float        phi()           const {this->okGP(); return phi_;} 
  float        r()             const {this->okGP(); return r_;}
  float        z()             const {this->okGP(); return z_;}
  float        rt()            const {this->okGP(); return rt_;}
  // %%% Those exclusively input to HT.
  float        dphi()          const {this->okHT(); this->valid(); return dphi_;}         // Delta(phi) = bend angle of stub
  float        rho()           const {this->okHT(); this->valid(); return rho_;}          // rho parameter
  float        phiS()          const {this->okHT(); return phiS_;}
  // Integer data after digitization (which doesn't degrade its resolution, but can recast it in a different form).
  // m bin range (= q/Pt bin range) allowed by bend filter
  // Note this range is centred on zero, so differs from Stub::min_qOverPt_bin() etc. which return a +ve number.
  int          m_min()         const {this->okHT(); return m_min_; }
  int          m_max()         const {this->okHT(); return m_max_; }
  // Tracker layer identifier encoded as it will be sent along optical link. 
  // Note that this differs from the encoding returned by Stub::layerIdReduced()!
  unsigned int iDigi_LayerID() const {this->okHT(); return iDigi_LayerID_;}
  // %%% Those exclusively input to GP.
  float        phiO()          const {this->okGP(); return phiO_;}
  float        bend()          const {this->okGP(); return bend_;}

  //--- The functions below give access to the original variables prior to digitization.
  //%%% Those common to GP & HT input.
  float        orig_phi()      const {this->okin(); return phi_orig_;} 
  float        orig_r()        const {this->okin(); return r_orig_;}
  float        orig_z()        const {this->okin(); return z_orig_;}
  //%%% Those exclusively input to HT.
  float        orig_dphi()     const {this->okin(); return dphi_orig_;}  // Delta(phi) = bend angle of stub
  float        orig_rho()      const {this->okin(); return rho_orig_;}   // rho parameter
  //%%% Those exclusively input to GP.
  float        orig_bend()     const {this->okin(); return bend_orig_;} 

  //--- Utility: return phi octant number corresponding to given phi sector number.
  unsigned int iGetOctant(unsigned int iPhiSec) const {return floor(iPhiSec*numPhiOctants_/numPhiSectors_);}


private:

  // Redigitize stub for input to Geographic Processor, if it was previously digitized for a different phi sector.
  void         quickMakeGPinput(int iPhiSec);

  // Redigitize stub for input to Hough transform, if it was previously digitized for a different phi sector.
  void         quickMakeHTinput(int iPhiSec);

  // Check that stub coords. & bend angle are within assumed digitization range.
  void         checkInRange() const;

  // Check that digitisation followed by undigitisation doesn't change significantly the stub coordinates.
  void         checkAccuracy() const;

  // Check that makeGPinput() or makeHTinput() are called before accessing digitized stub info.
  void         okGP()  const {if (! ranMakeGPinput_) throw cms::Exception("DigitalStub:: You forgot to call makeGPinput() or makeHTinput()!");}
  void         okHT()  const {if (! ranMakeHTinput_) throw cms::Exception("DigitalStub:: You forgot to call makeGPinput() or makeHTinput()!");}

  // Check that init() is called before accessing original pre-digitization variables.
  void         okin()  const {if (! ranInit_) throw cms::Exception("DigitalStub:: You forgot to call init()!");}

  // If using daisy-chain firmware, then it makes no sense to access the digiitzed values of dphi or rho.
  void         valid() const {if (iFirmwareType_ == 1) throw cms::Exception("DigitalStub:: You can't access digitized dphi or rho variables with daisy chain firmware!");}

private:

  //--- To check DigitialStub correctly initialized.
  bool                 ranInit_;
  bool                 ranMakeGPinput_;
  bool                 ranMakeHTinput_;

  //--- configuration

  // Digitization configuration
  int                  iFirmwareType_;
  unsigned int         phiSectorBits_;
  unsigned int         phiSBits_;
  float                phiSRange_;
  unsigned int         rtBits_;
  float                rtRange_;
  unsigned int         zBits_;
  float                zRange_;
  unsigned int         dPhiBits_;
  double               dPhiRange_;  // double precision owing to large number of digitisation bits.
  unsigned int         rhoBits_;
  float                rhoRange_;
  unsigned int         phiOBits_;
  double               phiORange_;
  unsigned int         bendBits_;
  float                bendRange_;

  // Digitization multipliers
  float                phiSMult_;
  float                rtMult_;
  float                zMult_;
  double               dPhiMult_; // double precision owing to large number of digitisation bits.
  float                rhoMult_;
  double               phiOMult_;
  float                bendMult_;

  // Are we using reduced layer ID, so layer can be packed into 3 bits?
  bool                 reduceLayerID_;

  // Number of phi sectors and phi octants.
  unsigned int         numPhiSectors_;
  unsigned int         numPhiOctants_;   // (=8)
  // Phi sector and phi octant width (radians)
  double               phiSectorWidth_;
  double               phiOctantWidth_;
  // Radius from beamline with respect to which stub r coord. is measured.
  float                chosenRofPhi_;

  // Number of q/Pt bins in Hough  transform array.
  unsigned int nbinsPt_;

  //--- Original floating point stub coords before digitization.
  float        phi_orig_;
  float        r_orig_;
  float        z_orig_;
  float        dphi_orig_;
  float        rho_orig_;
  float        rt_orig_;
  double       phiS_orig_;
  double       phiO_orig_;
  unsigned int layerID_; // Tracker layer ID
  unsigned int layerIDreduced_; // Tracker "reduced" layer ID
  unsigned int min_qOverPt_bin_orig_; // Range in q/Pt bins in HT array compatible with stub bend. (+ve definate)
  unsigned int max_qOverPt_bin_orig_; 
  float        bend_orig_;

  //--- Digits corresponding to stub coords.
  unsigned int iDigi_PhiSec_;
  int          iDigi_PhiS_;
  int          iDigi_Rt_;
  int          iDigi_Z_;
  int          iDigi_Dphi_;
  unsigned int iDigi_Rho_;
  unsigned int iDigi_LayerID_; // Encoded tracker layer
  int          m_min_; // Range in q/Pt bins in HT array compatible with stub bend. (range centred on zero)
  int          m_max_;
  unsigned int moduleType_;
  unsigned int iDigi_Octant_;
  int          iDigi_PhiO_;
  int          iDigi_Bend_;

  //--- Floating point stub coords derived from digitized info (so with degraded resolution).
  float phi_;
  float r_;
  float z_;
  float dphi_;
  float rho_;
  float phiS_;
  float rt_;
  float phiO_;
  float bend_;
};

#endif

