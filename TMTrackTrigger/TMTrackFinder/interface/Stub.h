#ifndef __STUB_H__
#define __STUB_H__

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"

// TTStubAssociationMap.h forgets to two needed files, so must include them here ...
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"

#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/DigitalStub.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include <vector>
#include <set>
#include <array>
#include <map>



class StackedTrackerGeometry;
class StackedTrackerDetId;
class TP;

typedef edm::Ref< edm::DetSetVector< PixelDigi >, PixelDigi > Ref_PixelDigi_;

typedef edmNew::DetSetVector< TTStub<Ref_PixelDigi_> > DetSetVec;
typedef edmNew::DetSet< TTStub<Ref_PixelDigi_> >       DetSet;
typedef edm::Ref<DetSetVec, TTStub<Ref_PixelDigi_> >   TTStubRef;
typedef edm::Ref< edmNew::DetSetVector< TTCluster< Ref_PixelDigi_ > >, TTCluster< Ref_PixelDigi_ > > TTClusterRef;
typedef TTStubAssociationMap<Ref_PixelDigi_>           TTStubAssMap;
typedef TTClusterAssociationMap<Ref_PixelDigi_>        TTClusterAssMap;



//=== Represents a Tracker stub (=pair of hits)

// TODO: remove me
// class Stub : public TTStubRef
class Stub
{
public:
	// Fill useful info about stub.
	Stub(TTStubRef ttStubRef, unsigned int index_in_vStubs, const Settings* settings, const StackedTrackerGeometry*  stackedGeometry);
	~Stub(){}

	bool operator==(const Stub& stubOther) {return (index() == stubOther.index());}

	// Fill truth info with association from stub to tracking particles.
	// The 1st argument is a map relating TrackingParticles to TP.
	void fillTruth(const std::map<edm::Ptr< TrackingParticle >, const TP* >& translateTP, edm::Handle<TTStubAssMap> mcTruthTTStubHandle, edm::Handle<TTClusterAssMap> mcTruthTTClusterHandle);

	// Calculate bin range along q/Pt axis of r-phi Hough transform array consistent with bend of this stub.
	void calcQoverPtrange();

	// Digitize stub for input to Geographic Processor, with digitized phi coord. measured relative to closest phi sector.
	// (This approximation is valid if their are an integer number of digitisation bins inside each phi octant).
	// However, you should also call digitizeForHTinput() before accessing digitized stub data, even if you only care about that going into GP! Otherwise, you will not identify stubs assigned to more than one octant.
	void digitizeForGPinput(unsigned int iPhiSec);

	// Digitize stub for input to Hough transform, with digitized phi coord. measured relative to specified phi sector.
	void digitizeForHTinput(unsigned int iPhiSec);

	// Restore stub to pre-digitized state. i.e. Undo what function digitize() did.

	void reset_digitize();

	// === Functions for returning info about reconstructed stubs ===

	// Location in InputData::vStubs_
	unsigned int                         index() const { return index_in_vStubs_; }

	//--- Stub data and quantities derived from it ---

	// Stub coordinates (optionally after digitisation, if digitisation requested via cfg).
	// N.B. Digitisation is not run when the stubs are created, but later, after stubs are assigned to sectors.
	// Until then, these functions return the original coordinates. 
	float                                  phi() const { return             phi_; }
	float                                    r() const { return               r_; }
	float                                    z() const { return               z_; }
	float                                 zTrk() const { return settings_->chosenRofZFilter()*z_/r_;}
	float                              zTrkRes() const { return  fabs(settings_->beamWindowZ()*(settings_->chosenRofZFilter() - r_)/r_) + fabs(settings_->chosenRofZFilter()*zErr_/r_) + fabs(settings_->chosenRofZFilter()*rErr_*z_/(r_*r_) );  }
	float                                  eta() const { return     asinh(z_/r_); }
	// Access to digitized version of stub coords.
	const DigitalStub&             digitalStub() const { return      digitalStub_;}
	// Access to booleans indicating if the stub has been digitized.
	bool                   digitizedForGPinput() const { return digitizedForGPinput_;}
	bool                   digitizedForHTinput() const { return digitizedForHTinput_;}

	// Get stub bend (i.e. displacement between two hits in stub in units of strip pitch) and its estimated resolution.
	float                                 bend() const { check1(); return bend_; } 
	// The bend resolution has a contribution from the sensor and a contribution from encoding the bend into
	// a reduced number of bits.
	float                              bendRes() const { return (settings_->bendResolution() + (numMergedBend_-1)*settings_->bendResolutionExtra()); }
	// Number of bend values which loss of bit to store bend resulted in being merged into this bend value.
	float                        numMergedBend() const { return numMergedBend_;}
	// Bend angle of track measured by stub and its estimated resolution.
	float                                 dphi() const { check2(); return dphi_; }
	float                              dphiRes() const { return (dphiOverBend() * bendRes()); }
	// Estimated track q/Pt based on stub bend info.
	float                              qOverPt() const { return (qOverPtOverBend() * bend()); }
	float                           qOverPtres() const { return (qOverPtOverBend() * bendRes()); }
	// Range in q/Pt bins in HT array compatible with stub bend.
	unsigned int               min_qOverPt_bin() const { return min_qOverPt_bin_; }
	unsigned int               max_qOverPt_bin() const { return max_qOverPt_bin_; }
	// Estimated phi0 of track at beam-line based on stub bend info.
	float                                 beta() const { return   (phi_ + dphi()); }
	// Estimated phi angle at which track crosses a given radius rad, based on stub bend info. Also estimate uncertainty on this angle due to endcap 2S module strip length. 
	// This is identical to beta() if rad=0.
	std::pair<float, float> trkPhiAtR(float rad) const;
	// Estimated resolution in trkPhiAtR(rad) based on nominal stub bend resolution.
	float                trkPhiAtRres(float rad) const { return dphiRes() * fabs(1 - rad / r_); }
	// Difference in phi between stub and angle at which track crosses given radius, assuming track has given Pt.
	float           phiDiff(float rad, float Pt) const { return fabs(r_ - rad)*(settings_->invPtToDphi())/Pt; }
	// -- conversion factors
	// Ratio of bend angle to bend, where bend is the displacement in strips between the two hits making up stub.
	float                         dphiOverBend() const { check2(); return dphiOverBend_; }
	// Two related parameters defined in spec. document, one of which is passed from PP to MP along optical link.
	float                      rhoRawParameter() const { return dphiOverBend();} // A pseudonym of dPhiOverBend.
	float                         rhoParameter() const { return rhoRawParameter() * bendRes();} 
	// Ratio of q/Pt to bend, where bend is the displacement in strips between the two hits making up stub.
	float                      qOverPtOverBend() const { return dphiOverBend() / (r_ * settings_->invPtToDphi()); }

	//--- Info about the two clusters that make up the stub.
	// Coordinates in frame of sensor, measured in units of strip pitch along two orthogonal axes running perpendicular and parallel to longer axis of pixels/strips (U & V).
	std::array<float, 2>        localU_cluster() const { return localU_cluster_;}
	std::array<float, 2>        localV_cluster() const { return localV_cluster_;}

	//--- Check if this stub will be output by front-end readout electronics,
	//--- (where we can reconfigure the stub window size and rapidity cut).
	//--- Don't use stubs failing this cut.
	bool                          frontendPass() const { return    frontendPass_; }
	// Indicates if stub would have passed front-end cuts, were it not for window size encoded in DataCorrection.h
	bool              stubFailedDataCorrWindow() const { return    stubFailedDataCorrWindow_;}

	//--- Quantities common to all stubs in a given module ---

	// Unique identifier for each stacked module, allowing one to check which stubs are on the same module.
	unsigned int                         idDet() const { return           idDet_;}
	// Uncertainty in stub coordinates due to strip length, assumed equal to 0.5*strip-length in 2S modules and zero in PS modules.
	float                                 rErr() const { return            rErr_;}
	float                                 zErr() const { return            zErr_;}
	// Coordinates of centre of two sensors in (r,phi,z)
	float                                 minR() const { return      moduleMinR_; }
	float                                 maxR() const { return      moduleMaxR_; }
	float                               minPhi() const { return    moduleMinPhi_; }
	float                               maxPhi() const { return    moduleMaxPhi_; }
	float                                 minZ() const { return      moduleMinZ_; }
	float                                 maxZ() const { return      moduleMaxZ_; }
	// Sensor pitch over separation.
	float                         pitchOverSep() const { return    pitchOverSep_;}
	// Location of stub in module in units of strip number (or pixel number along finest granularity axis).
	// Range from 0 to (nStrips - 1) inclusive.
	unsigned int                          iphi() const { return            iphi_; }
	// Module type: PS or 2S?
	bool                              psModule() const { return        psModule_;}
	// Tracker layer ID number (1-6 = barrel layer; 11-15 = endcap A disk; 21-25 = endcap B disk)
	unsigned int                       layerId() const { return         layerId_; }
	// Reduced layer ID (in range 1-7). This encodes the layer ID in only 3 bits (to simplify firmware) by merging some barrel layer and endcap disk layer IDs into a single ID.
	unsigned int                layerIdReduced() const;
	// Endcap ring of module (returns zero in case of barrel)
	unsigned int                    endcapRing() const { return      endcapRing_; }
	bool                                barrel() const { return          barrel_; }

	// Strip pitch (or pixel pitch along shortest axis).
	float                           stripPitch() const { return      stripPitch_; } 
	// Strip length (or pixel pitch along longest axis).
	float                          stripLength() const { return     stripLength_; } 
	// No. of strips in sensor.
	unsigned int                       nStrips() const { return         nStrips_; }
	// Width of sensitive region of sensor.
	float                          sensorWidth() const { return     sensorWidth_; }
	// Hit resolution perpendicular to strip (or to longest pixel axis) = pitch/sqrt(12). Measures phi.
	float                            sigmaPerp() const { return       sigmaPerp_; }
	// Hit resolution parallel to strip (or to longest pixel axis) = length/sqrt(12). Measures r or z.
	float                             sigmaPar() const { return        sigmaPar_; }

	// Clone a few of the above functions with the less helpful names expected by the track fitting code. (Try to phase these out with time ...)
	unsigned int                        nstrip() const { return nStrips(); }
	float                                width() const { return sensorWidth(); }
	float                               sigmaX() const { return sigmaPerp(); }
	float                               sigmaZ() const { return sigmaPar(); }

	//--- Truth info

	// Association of stub to tracking particles
	std::set<const TP*>               assocTPs() const { return        assocTPs_; } // Return TPs associated to this stub. (Whether only TPs contributing to both clusters are returned is determined by "StubMatchStrict" config param.)
	bool	 			     genuine() const { return (assocTPs_.size() > 0); } // Did stub match at least one TP?
	const TP*                          assocTP() const { return         assocTP_; } // If only one TP contributed to both clusters, this tells you which TP it is. Returns nullptr if none.

	// Association of both clusters making up stub to tracking particles
	std::array<bool, 2>	        genuineCluster() const { return std::array<bool, 2>{ {(assocTPofCluster_[0] != nullptr), (assocTPofCluster_[1] != nullptr)} }; } // Was cluster produced by a single TP?
	std::array<const TP*, 2>  assocTPofCluster() const { return       assocTPofCluster_; } // Which TP made each cluster. Warning: If cluster was not produced by a single TP, then returns nullptr! (P.S. If both clusters match same TP, then this will equal assocTP()).

	// Note if stub is a crazy distance from the tracking particle trajectory that produced it. (e.g. perhaps produced by delta ray)
	bool                             crazyStub() const;

	// Get stub bend and its resolution, as available within the front end chip (i.e. prior to loss of bits
	// or digitisation).
	float                       bendInFrontend() const { return bendInFrontend_; } 
	float                    bendResInFrontend() const { return settings_->bendResolution(); } 

	TTStubRef                   cmsswTTStubRef() const { return cmssswTTStubRef_; }
	
private:

	// Degrade assumed stub bend resolution.
	// Also return boolean indicating if stub bend was outside assumed window, so stub should be rejected
	// and return an integer indicating how many values of bend are merged into this single one.
	void degradeResolution(float bend, const StackedTrackerDetId& stDetId,
						float& degradedBend, bool& reject, unsigned int& num) const;

	// Set the frontendPass_ flag, indicating if frontend readout electronics will output this stub.  
	// Argument indicates if stub bend was outside window size encoded in DataCorrection.h
	void setFrontend(bool rejectStub);          

	// Set info about the module that this stub is in.
	void  setModuleInfo(const StackedTrackerGeometry* stackedGeometry, const StackedTrackerDetId& stDetId);

	// Function to set rho parameter value. Since the rho parameter is not a data member of this class, this is done by setting the value of 
	// dphiOverBend_ (also known as "raw rho") from which rho is derived.
	void  setRhoParameter(float rho) { dphiOverBend_ = rho / bendRes(); }

	// No HT firmware can access directly the stub bend info.
	void  check1() const {if (digitizedForHTinput_) throw cms::Exception("DigitalStub:: You can't access digitized bend variable within HT firmware!");}
	// If using daisy-chain firmware, then it makes no sense to access the digiitzed values of dphi or rho within HT.
	void  check2() const {if (digitizedForHTinput_ && settings_->firmwareType() == 1) throw cms::Exception("DigitalStub:: You can't access digitized dphi or rho variables within daisy chain HT firmware!");}

private:

	const Settings* settings_; // configuration parameters.

	unsigned int                     index_in_vStubs_; // location of this stub in InputData::vStubs

	//--- Parameters passed along optical links from PP to MP (or equivalent ones if easier for analysis software to use).
	// N.B. Parameters dphiOverBend_ and dphi_ are used with the systolic & 2-c-bin firmware, whilst parameters
	// min_qOverPt_bin_ & max_qOverPt_bin_ are used with the daisy-chain firmware.
	// WARNING: If you add any variables in this section, take care to ensure that they are digitized correctly by Stub::digitize().
	float                                        phi_; // stub coords, optionally after digitisation.
	float                                          r_;
	float                                          z_;
	float                                       bend_; // bend of stub.
	float                               dphiOverBend_; // related to rho parameter.
	float                                       dphi_;
	unsigned int                     min_qOverPt_bin_; // Range in q/Pt bins in HT array compatible with stub bend.
	unsigned int                     max_qOverPt_bin_; 

	//--- Info about the two clusters that make up the stub.
	std::array<float, 2>              localU_cluster_;
	std::array<float, 2>              localV_cluster_;

	//--- Parameters common to all stubs in a given module.
	unsigned int                               idDet_; 
	float                                       rErr_;
	float                                       zErr_;
	float                                 moduleMinR_;
	float                                 moduleMaxR_;
	float                               moduleMinPhi_;
	float                               moduleMaxPhi_;
	float                                 moduleMinZ_;
	float                                 moduleMaxZ_;
	float                               pitchOverSep_;
	unsigned int                                iphi_;
	bool                                    psModule_;
	unsigned int                             layerId_;
	unsigned int                          endcapRing_;
	bool                                      barrel_;
	float                                  sigmaPerp_;
	float                                   sigmaPar_;
	float                                 stripPitch_;
	float                                stripLength_;
	unsigned int                             nStrips_;
	float                                sensorWidth_;

	//--- Truth info about stub.
	const TP*                                assocTP_;
	std::set<const TP*>                     assocTPs_;
	//--- Truth info about the two clusters that make up the stub
	std::array<const TP*, 2>        assocTPofCluster_;

	// Would front-end electronics output this stub?
	bool                                frontendPass_;
	// Did stub fail window cuts assumed in DataCorrection.h?
	bool                    stubFailedDataCorrWindow_;
	// Bend in front end chip (prior to degredation by loss of bits & digitization).
	float                             bendInFrontend_;
	// Used for stub bend resolution degrading.
	unsigned int                       numMergedBend_;

	DigitalStub                          digitalStub_; // Class used to digitize stub if required.
	bool                         digitizedForGPinput_; // Has this stub been digitized for GP input?
	bool                         digitizedForHTinput_; // Has this stub been digitized for HT input?
	
	TTStubRef                        cmssswTTStubRef_;
};

#endif
