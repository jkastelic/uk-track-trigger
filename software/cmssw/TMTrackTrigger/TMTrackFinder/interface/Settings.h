#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <vector>
#include <iostream>

using namespace std;

// Stores all configuration parameters + some hard-wired constants.

class Settings {

public:
  Settings(const edm::ParameterSet& iConfig);
  ~Settings(){}

  //=== Cuts on MC truth tracks for tracking efficiency measurements.

  double               genMinPt()                const   {return genMinPt_;}
  double               genMaxAbsEta()            const   {return genMaxAbsEta_;}
  double               genMaxVertR()             const   {return genMaxVertR_;}
  double               genMaxVertZ()             const   {return genMaxVertZ_;}
  vector<int>          genPdgIds()               const   {return genPdgIds_;}
  // Additional cut on MC truth tracks for algorithmic tracking efficiency measurements.
  unsigned int         genMinStubLayers()        const   {return genMinStubLayers_;} // Min. number of layers TP made stub in.

  //=== Cuts applied to stubs before arriving in L1 track finding board.

  // Reduce number of bits used by front-end chips to store stub bend info.
  bool                 bendResReduced()          const   {return bendResReduced_;}
  // Don't use stubs with eta beyond this cut, since the tracker geometry makes it impossible to reconstruct tracks with them.
  double               maxStubEta()              const   {return maxStubEta_;}
  // Don't use stubs whose measured Pt from bend info is significantly below HTArraySpec.HoughMinPt, where "significantly" means allowing for resolution in q/Pt derived from stub bend resolution HTFilling.BendResolution
  bool                 killLowPtStubs()          const   {return killLowPtStubs_;}
  // Bend resolution assumed by bend filter in units of strip pitch. Also used when assigning stubs to sectors if calcPhiTrkRes() is true.
  double               bendResolution()          const   {return bendResolution_;}        
  // Additional contribution to bend resolution from its encoding into a reduced number of bits.
  // This number is the assumed resolution relative to the naive guess of its value.
  // It is ignored in BendResReduced = False.
  double               bendResolutionExtra()     const   {return bendResolutionExtra_;}        

  //=== Optional stub digitization configuration

  bool                 enableDigitize()          const   {return enableDigitize_;}
  unsigned int         firmwareType()            const   {return firmwareType_;}
  //--- Parameters available in MP board.
  unsigned int         phiSectorBits()           const   {return phiSectorBits_;}
  unsigned int         phiSBits()                const   {return phiSBits_;}
  double               phiSRange()               const   {return phiSRange_;}
  unsigned int         rtBits()                  const   {return rtBits_;}
  double               rtRange()                 const   {return rtRange_;}
  unsigned int         zBits()                   const   {return zBits_;}
  double               zRange()                  const   {return zRange_;}
  unsigned int         dPhiBits()                const   {return dPhiBits_;}
  double               dPhiRange()               const   {return dPhiRange_;}
  unsigned int         rhoBits()                 const   {return rhoBits_;}
  double               rhoRange()                const   {return rhoRange_;}
  //--- Parameters available in GP board (excluding any in common with MP specified above).
  unsigned int         phiOBits()                const   {return phiOBits_;}
  double               phiORange()               const   {return phiORange_;}
  unsigned int         bendBits()                const   {return bendBits_;}

  //=== Definition of phi sectors.

  unsigned int         numPhiSectors()           const   {return numPhiSectors_;}
  double               chosenRofPhi()            const   {return chosenRofPhi_;} // Use phi of track at this radius for assignment of stubs to phi sectors & also for one of the axes of the r-phi HT. If ChosenRofPhi=0, then use track phi0.
  bool                 useStubPhi()              const   {return useStubPhi_;} // Require stub phi to be consistent with track of Pt > HTArraySpec.HoughMinPt that crosses HT phi axis?
  bool                 useStubPhiTrk()           const   {return useStubPhiTrk_;} // Require stub phi0 (or phi65 etc.) as estimated from stub bend, to lie within HT phi axis, allowing tolerance specified below?
  double               assumedPhiTrkRes()        const   {return assumedPhiTrkRes_;} // Tolerance in stub phi0 (or phi65) assumed to be this fraction of phi sector width. (N.B. If > 0.5, then stubs can be shared by more than 2 phi sectors).
  bool                 calcPhiTrkRes()           const   {return calcPhiTrkRes_;} // If true, tolerance in stub phi0 (or phi65 etc.) will be reduced below AssumedPhiTrkRes if stub bend resolution specified in StubCuts.BendResolution suggests it is safe to do so.
  bool                 handleStripsPhiSec()      const   {return handleStripsPhiSec_;} // Should algorithm allow for uncertainty in stub (r,z) coordinate caused by length of 2S module strips when assigning stubs to phi sectors?

  //=== Definition of eta sectors.

  vector<double>       etaRegions()              const   {return etaRegions_;} // Boundaries of eta regions de
  unsigned int         numEtaRegions()           const   {return (etaRegions_.size() - 1);}
  double               chosenRofZ()              const   {return chosenRofZ_;} // Use z of track at this radius for assignment of stubs to phi sectors & also for one of the axes of the r-z HT. 
  double               beamWindowZ()             const   {return beamWindowZ_;} // Half-width of window supposed to contain beam-spot in z.
  bool                 handleStripsEtaSec()      const   {return handleStripsEtaSec_;} // Should algorithm allow for uncertainty in stub (r,z) coordinate caused by length of 2S module strips when assigning stubs to eta sectors?

  //=== r-phi Hough transform array specifications.
                                
  double               houghMinPt()              const   {return houghMinPt_;}
  unsigned int         houghNbinsPt()            const   {return houghNbinsPt_;} // Dimension in any q/Pt related variable. Not valid if houghNcellsRphi() > 0.
  unsigned int         houghNbinsPhi()           const   {return houghNbinsPhi_;} // Dimension in any track-phi related variable. Not valid if houghNcellsRphi() > 0.
  int                  houghNcellsRphi()         const   {return houghNcellsRphi_;} // Required no. of cells in r-phi HT array. If > 0, then parameters HoughNbinsPt and HoughNbinsPhi will be calculated from the constraints that their product should equal HoughNcellsRz and their ratio should make the maximum |gradient|" of stub lines in the HT array equal to 1. If <= 0, then HoughNbinsPt and HoughNbinsPhi will be taken from the values configured above. 
  bool                 enableMerge2x2()          const   {return enableMerge2x2_;} // Groups of neighbouring 2x2 cells in HT will be treated as if they are a single large cell.
  double               maxPtToMerge2x2()         const   {return maxPtToMerge2x2_;} // but only cells with pt < maxPtToMerge2x2() will be merged in this way (irrelevant if enableMerge2x2() = false).
  unsigned int         numSubSecsEta()           const   {return numSubSecsEta_;} // Subdivide each sector into this number of subsectors in eta within r-phi HT.

  //=== r-z Hough transform array specifications.
                                
  bool                 enableRzHT()              const   {return enableRzHT_;} // Use the r-z Hough Transform?
  unsigned int         houghNbinsZ0()            const   {return houghNbinsZ0_;} // Dimension in z0. Not valid if houghNcellsRz() > 0.
  unsigned int         houghNbinsZ65()           const   {return houghNbinsZ65_;} // Dimension in z65 (or any other z related variable). Not valid if houghNcellsRz() > 0.
  int                  houghNcellsRz()           const   {return houghNcellsRz_;} // Required no. of cells in r-z HT array. If > 0, then parameters HoughNbinsZ0 and HoughNbinsZ65 will be calculated from the constraints that their product should equal HoughNcellsRz and their ratio should make the maximum |gradient|" of stub lines in the HT array equal to 1. If <= 0, then HoughNbinsZ0 and HoughNbinsRz will be taken from the values configured above. 

  //=== Rules governing how stubs are filled into the r-phi Hough Transform array.
                                
  // Should algorithm allow for uncertainty in stub (r,z) coordinate caused by length of 2S module strips when filling stubs in r-phi HT?
  bool                 handleStripsRphiHT()      const   {return handleStripsRphiHT_;} 
  // Take all cells in HT array crossed by line corresponding to each stub (= 0) or take only some to reduce rate at cost 
  // of efficiency ( > 0). If this option is > 0, it can be 1 or 2, corresponding to different algorithms for rejecting some of the cells. 
  unsigned int         killSomeHTCellsRphi()     const   {return killSomeHTCellsRphi_;} 
  // Use filter in each HT cell using only stubs which have consistent bend, allowing for resolution specified in StubCuts.BendResolution.
  bool                 useBendFilter()           const   {return useBendFilter_;} 
  // A filter is used each HT cell, which prevents more than the specified number of stubs being stored in the cell. (Reflecting memory limit of hardware).   
  unsigned int         maxStubsInCell()          const   {return maxStubsInCell_;}
  // If this returns true, and if more than busySectorNumStubs() stubs are assigned to tracks by an r-phi HT array, then 
  // the excess tracks are killed, with lowest Pt ones killed first. This is because hardware has finite readout time.
  bool                 busySectorKill()          const   {return busySectorKill_;}
  unsigned int         busySectorNumStubs()      const   {return busySectorNumStubs_;}
  // If this is True, then the BusySectorNumStubs cut is applied to +ve and -ve charge track seperately. (Irrelevant if BusySectorKill = False).
  bool                 busySectorEachCharge()    const   {return busySectorEachCharge_;}

  //=== Rules governing how stubs are filled into the r-z Hough Transform array. (Irrelevant if enableRzHT = false.)
                                
  // Should algorithm allow for uncertainty in stub (r,z) coordinate caused by length of 2S module strips when filling stubs in r-z HT?
  bool                 handleStripsRzHT()        const   {return handleStripsRzHT_;} 
  // Take all cells in HT array crossed by line corresponding to each stub (= 0) or take only some to reduce rate at cost 
  // of efficiency ( > 0). If this option is > 0, it can be 1 or 2, corresponding to different algorithms for rejecting some of the cells. 
  unsigned int         killSomeHTCellsRz()       const   {return killSomeHTCellsRz_;} 

  //=== Options controlling r-z track filters (or any other track filters run after the Hough transform, as opposed to inside it).

  // Use filter in each HT cell using only stubs which have consistent rapidity?
  bool                 useEtaFilter()            const   {return useEtaFilter_;}
  // Use filter in each HT cell using only stubs which have consistent value of ZTrk?
  bool                 useZTrkFilter()           const   {return useZTrkFilter_;}
  // Filter Stubs in each HT cell using a tracklet seed algorithm?
  bool                 useSeedFilter()           const   {return useSeedFilter_;}
  // Use z of track at this radius (relevant for Ztrk filter). 
  double               chosenRofZFilter()        const   {return chosenRofZFilter_;}
  // Added resolution beyond that estimated from hit resolution (relevant for Seed filter).
  double               seedResolution()          const   {return seedResolution_;}     
  // Store stubs compatible with all possible good seed (relevant for Seed filter)?
  bool                 keepAllSeed()             const   {return keepAllSeed_;}
  // Maximum number of seed combinations to check (relevant for Seed filter).
  unsigned int         maxSeedCombinations()     const   {return maxSeedCombinations_;}
  // Check that estimated zTrk from seeding stub is within the sector boundaries (relevant for Seed filter)?
  bool                 zTrkSectorCheck()         const   {return zTrkSectorCheck_;}       

  //=== Rules for deciding when the track finding has found an L1 track candidate

  // Min. number of layers in HT cell that must have stubs for track to be declared found.
  unsigned int         minStubLayers()           const   {return minStubLayers_;}
  // Reduce MinStubLayers by 1 for HT cells corresponding to track pt above this cut.
  // If this is set to >= 10000., this option is disabled.
  double               minPtToReduceLayers()     const   {return minPtToReduceLayers_;}
  // Define layers using layer ID (true) or by bins in radius of 5 cm width (false)?
  bool                 useLayerID()              const   {return useLayerID_;}
  //Reduce this layer ID, so that it takes no more than 8 different values in any eta region (simplifies firmware)?
  bool                 reduceLayerID()           const   {return reduceLayerID_;}

  //=== Specification of algorithm to eliminate duplicate tracks

  // --- Which algorithms are enabled.
  // Algorithm used for duplicate removal of 2D tracks produced by r-phi HT.
  unsigned int         dupTrkAlgRphi()           const   {return dupTrkAlgRphi_;}
  // Algorithm used for duplicate removal run on 2D tracks produced by r-z HT corresponding to a single r-phi HT track.
  unsigned int         dupTrkAlgRz()             const   {return dupTrkAlgRz_;}
  // Algorithm run on all 3D tracks within each sector after r-z HT or r-z seed filter.
  unsigned int         dupTrkAlgRzSeg()          const   {return dupTrkAlgRzSeg_;}
  // Algorithm run on tracks after the track helix fit has been done.
  unsigned int         dupTrkAlgFit()            const   {return dupTrkAlgFit_;}
  //--- Options used by individual algorithms.
  unsigned int         dupTrkMinIndependent()    const   {return dupTrkMinIndependent_;}
  unsigned int         dupTrkMinCommonHitsLayers() const {return dupTrkMinCommonHitsLayers_;}
  double               dupTrkChiSqCut()           const  {return dupTrkChiSqCut_;}
  double               dupMaxQOverPtScan()       const   {return dupMaxQOverPtScan_;}
  double               dupMaxPhi0Scan()          const   {return dupMaxPhi0Scan_;}
  double               dupMaxZ0Scan()            const   {return dupMaxZ0Scan_;}
  double               dupMaxTanLambdaScan()     const   {return dupMaxTanLambdaScan_;}

  //=== Rules for deciding when a reconstructed L1 track matches a MC truth particle (i.e. tracking particle).

  //--- Three different ways to define if a tracking particle matches a reco track candidate. (Usually, set two of them to ultra loose).
  // Min. fraction of matched stubs relative to number of stubs on reco track.
  double               minFracMatchStubsOnReco() const   {return minFracMatchStubsOnReco_;}
  // Min. fraction of matched stubs relative to number of stubs on tracking particle.
  double               minFracMatchStubsOnTP()   const   {return minFracMatchStubsOnTP_;}
  // Min. number of matched layers.
  unsigned int         minNumMatchLayers()       const   {return minNumMatchLayers_;}
  // Associate stub to TP only if the TP contributed to both its clusters? (If False, then associate even if only one cluster was made by TP).
  bool                 stubMatchStrict()         const   {return stubMatchStrict_;}

  //=== Track Fitting Settings

  //--- Options applicable to all track fitters ---

  // Track fitting algorithms to use. You can run several in parallel.
  std::vector<std::string> trackFitters()        const   {return trackFitters_;}
  // Cut on chi2/dof of fitted track.
  double               chi2OverNdfCut()          const   {return chi2OverNdfCut_;}
  // Print detailed summary of track fit performance at end of job (as opposed to a brief one)?
  bool                 detailedFitOutput()       const   {return detailedFitOutput_;} 

  //--- Options for chi2 & Linear Regression track fitters ---

  // Number of iterations that the track fit should perform.
  unsigned int         numTrackFitIterations()   const   {return numTrackFitIterations_;}
  // Optionally remove hit with worst residual in track fit? (Only used by chi2 track fit).
  bool                 killTrackFitWorstHit()    const   {return killTrackFitWorstHit_;}
  // Cuts in standard deviations used to kill hits with big residuals during fit. If the residual exceeds the "General" 
  // cut, the hit is killed providing it leaves the track with enough hits to survive. If the residual exceeds the 
  // "Killing" cut, the hit is killed even if that kills the track.
  double               generalResidualCut()      const   {return generalResidualCut_;}
  double               killingResidualCut()      const   {return killingResidualCut_;}

  //--- Additional options for Linear Regression track fitter ---

  bool                 combineResiduals()               const { return combineResiduals_; }
  bool                 lineariseStubPosition()          const { return lineariseStubPosition_; }
  bool                 checkSectorConsistency()         const { return checkSectorConsistency_; }
  bool                 checkHTCellConsistency()         const { return checkHTCellConsistency_; }
  unsigned int         minPSLayers()                    const { return minPSLayers_; }
  bool                 fitPerfectCandidatesOnly()       const { return fitPerfectCandidatesOnly_; }

  //--- Options for Kalman filter track fitters ---

  // larger number has more debugging outputs.
  unsigned             kalmanDebugLevel()               const { return kalmanDebugLevel_; }
  // Internal histograms are filled if it is True
  bool                 kalmanFillInternalHists()        const { return kalmanFillInternalHists_;}
  // Multiple scattering factor.  Not working. Set to 0.
  double               kalmanMultiScattFactor()         const { return kalmanMultiScattFactor_;}
  // A stub which is inconsistent with the state is not processed for the state. Cut value on chisquare from the forcast is set.
  double               kalmanValidationGateCutValue()   const { return kalmanValidationGateCutValue_; }
  // Best candidate is selected from the candidates with the more number of stubs if this is True. Chi2 only selection is better to remove combinatorial hits
  bool                 kalmanSelectMostNumStubState()   const { return kalmanSelectMostNumStubState_; }
  // After this number of stubs in the state, the rest of the stubs are not processed and go to the next state.
  unsigned             kalmanMaxNumNextStubs()          const { return kalmanMaxNumNextStubs_; }
  // Allowed number of virtual stubs. (CAUTION: kalmanState::nVirtualStubs() return the # of virtual stubs including the seed. The number on this option does not include the seed.)
  unsigned             kalmanMaxNumVirtualStubs()       const { return kalmanMaxNumVirtualStubs_; }
  // The number of states for each virtual stub list is restricted to this number.  The lists are truncated to this number of states.
  unsigned             kalmanMaxNumStatesCutValue()     const { return kalmanMaxNumStatesCutValue_; } 
  // The state is removed if the reduced chisquare is more than this number.
  double               kalmanStateReducedChi2CutValue() const { return kalmanStateReducedChi2CutValue_; }

  //=== Debug printout
  unsigned int         debug()                   const   {return debug_;}

  // Booleain indicating if an output EDM file will be written.
  // N.B. This parameter does not appear inside TMTrackProducer_Defaults_cfi.py . It is created inside tmtt_tf_analysis_cfg.py .
  bool                 writeOutEdmFile()         const   {return writeOutEdmFile_;}

  //=== Hard-wired constants
  double               pitchPS()                 const   {cout<<"ERROR: Use Stub::stripPitch instead of Settings::pitchPS!";exit(1);return 0.;} // pitch of PS modules - OBSOLETE
  double               pitch2S()                 const   {cout<<"ERROR: Use Stub::stripPitch instead of Settings::pitch2S!";exit(1);return 0.;} // pitch of 2S modules - OBSOLETE
  double               invPtToInvR()             const   {return (this->getBfield())*(3.0E8/1.0E11);} // B*c/1E11 - converts q/Pt to 1/radius_of_curvature
  double               invPtToDphi()             const   {return (this->getBfield())*(3.0E8/2.0E11);} // B*c/2E11 - converts q/Pt to track angle at some radius from beamline.  
  double               trackerOuterRadius()      const   {return 109.2;}  // max. occuring stub radius.
  double               trackerInnerRadius()      const   {return  21.3;}  // min. occuring stub radius.
  double               trackerHalfLength()       const   {return 270.;}  // half-length of tracker. 
  double               stripLength2S()           const   {cout<<"ERROR: Use Stub::stripLength instead of Settings::stripLength2S!"<<endl;exit(1);return 0.;}    // Strip length of 2S modules. - OBSOLETE
  double               layerIDfromRadiusBin()    const   {return 6.;}    // When counting stubs in layers, actually histogram stubs in distance from beam-line with this bin size.
  double               crazyStubCut()            const   {return 0.01;}  // Stubs differing from TP trajectory by more than this in phi are assumed to come from delta rays etc.

  //=== Set and get B-field value in Tesla. 
  // N.B. This must bet set for each event, and can't be initialized at the beginning of the job.
  void                 setBfield(float bField)           {bField_ = bField;}
  float                getBfield()               const   {if (bField_ == 0.) throw cms::Exception("Settings.h:You attempted to access the B field before it was initialized"); return bField_;}

private:

  // Parameter sets for differents types of configuration parameter.
  edm::ParameterSet    genCuts_;
  edm::ParameterSet    stubCuts_;
  edm::ParameterSet    stubDigitize_;
  edm::ParameterSet    phiSectors_;
  edm::ParameterSet    etaSectors_;
  edm::ParameterSet    htArraySpecRphi_;
  edm::ParameterSet    htArraySpecRz_;
  edm::ParameterSet    htFillingRphi_;
  edm::ParameterSet    htFillingRz_;
  edm::ParameterSet    rzFilterOpts_;
  edm::ParameterSet    l1TrackDef_;
  edm::ParameterSet    dupTrkRemoval_;
  edm::ParameterSet    trackMatchDef_;
  edm::ParameterSet    trackFitSettings_;

  // Cuts on truth tracking particles.
  double               genMinPt_;
  double               genMaxAbsEta_;
  double               genMaxVertR_;
  double               genMaxVertZ_;
  vector<int>          genPdgIds_;
  unsigned int         genMinStubLayers_;

  // Cuts applied to stubs before arriving in L1 track finding board.
  bool                 bendResReduced_;
  double               maxStubEta_;
  bool                 killLowPtStubs_;
  double               bendResolution_;
  double               bendResolutionExtra_;

  // Optional stub digitization.
  bool                 enableDigitize_;
  unsigned int         firmwareType_;
  unsigned int         phiSectorBits_;
  unsigned int         phiSBits_;
  double               phiSRange_;
  unsigned int         rtBits_;
  double               rtRange_;
  unsigned int         zBits_;
  double               zRange_;
  unsigned int         dPhiBits_;
  double               dPhiRange_;
  unsigned int         rhoBits_;
  double               rhoRange_;
  unsigned int         phiOBits_;
  double               phiORange_;
  unsigned int         bendBits_;

  // Definition of phi sectors.
  unsigned int         numPhiSectors_;
  double               chosenRofPhi_;  
  bool                 useStubPhi_;
  bool                 useStubPhiTrk_;
  double               assumedPhiTrkRes_;
  bool                 calcPhiTrkRes_;
  bool                 handleStripsPhiSec_;

  // Definition of eta sectors.
  vector<double>       etaRegions_;                                
  double               chosenRofZ_;  
  double               beamWindowZ_;
  bool                 handleStripsEtaSec_;
                                
  // r-phi Hough transform array specifications.
  double               houghMinPt_;
  unsigned int         houghNbinsPt_;
  unsigned int         houghNbinsPhi_;
  int                  houghNcellsRphi_;
  bool                 enableMerge2x2_;
  double               maxPtToMerge2x2_;
  unsigned int         numSubSecsEta_;

  // r-z Hough transform array specifications.
  bool                 enableRzHT_;
  unsigned int         houghNbinsZ0_;
  unsigned int         houghNbinsZ65_;
  int                  houghNcellsRz_;
                                
  // Rules governing how stubs are filled into the r-phi Hough Transform array.
  bool                 handleStripsRphiHT_;
  unsigned int         killSomeHTCellsRphi_;
  bool                 useBendFilter_;
  unsigned int         maxStubsInCell_;
  bool                 busySectorKill_;
  unsigned int         busySectorNumStubs_;
  bool                 busySectorEachCharge_; 
  
  // Rules governing how stubs are filled into the r-z Hough Transform array. (Irrelevant if enableRzHT = false.)
  bool                 handleStripsRzHT_;
  unsigned int         killSomeHTCellsRz_;

  // Options controlling r-z track filters (or any other track filters run after the Hough transform, as opposed to inside it).
  bool                 useEtaFilter_;
  bool                 useZTrkFilter_;
  bool                 useSeedFilter_;
  double               chosenRofZFilter_;
  double               seedResolution_;
  bool                 keepAllSeed_;
  unsigned int         maxSeedCombinations_;
  bool                 zTrkSectorCheck_;

  // Rules for deciding when the track-finding has found an L1 track candidate
  unsigned int         minStubLayers_;
  double               minPtToReduceLayers_;
  bool                 useLayerID_;
  bool                 reduceLayerID_;

  // Specification of algorithm to eliminate duplicate tracks
  unsigned int         dupTrkAlgRphi_;
  unsigned int         dupTrkAlgRz_;
  unsigned int         dupTrkAlgRzSeg_;
  unsigned int         dupTrkAlgFit_;
  unsigned int         dupTrkMinIndependent_;
  unsigned int         dupTrkMinCommonHitsLayers_;
  double               dupTrkChiSqCut_;
  double               dupMaxQOverPtScan_;
  double               dupMaxPhi0Scan_;
  double               dupMaxZ0Scan_;
  double               dupMaxTanLambdaScan_;

  // Rules for deciding when a reconstructed L1 track matches a MC truth particle (i.e. tracking particle).
  double               minFracMatchStubsOnReco_;
  double               minFracMatchStubsOnTP_;
  unsigned int         minNumMatchLayers_;
  bool                 stubMatchStrict_;

  // Track Fitting Settings
  std::vector<std::string> trackFitters_;
  double               chi2OverNdfCut_;
  bool                 detailedFitOutput_;
  unsigned int         numTrackFitIterations_;
  bool                 killTrackFitWorstHit_;
  double               generalResidualCut_;
  double               killingResidualCut_;
  bool                 combineResiduals_;
  bool                 lineariseStubPosition_;
  bool                 checkSectorConsistency_;
  bool                 checkHTCellConsistency_;
  unsigned int         minPSLayers_;
  bool                 fitPerfectCandidatesOnly_;
  unsigned             kalmanDebugLevel_;
  bool                 kalmanFillInternalHists_;
  double               kalmanMultiScattFactor_; 
  double               kalmanValidationGateCutValue_;
  bool                 kalmanSelectMostNumStubState_;
  unsigned             kalmanMaxNumNextStubs_;
  unsigned             kalmanMaxNumVirtualStubs_;
  unsigned             kalmanMaxNumStatesCutValue_;
  double               kalmanStateReducedChi2CutValue_;
  
  // Debug printout
  unsigned int         debug_;

  // Boolean indicating an an EDM output file will be written.
  bool                 writeOutEdmFile_;

  // B-field in Tesla
  float                bField_;
};

#endif
