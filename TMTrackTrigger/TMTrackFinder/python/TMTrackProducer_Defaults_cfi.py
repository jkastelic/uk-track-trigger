import FWCore.ParameterSet.Config as cms

TMTrackProducer = cms.EDProducer('TMTrackProducer',

  #=== Cuts on MC truth particles (i.e., tracking particles) used for tracking efficiency measurements.

  GenCuts = cms.PSet(
     GenMinPt         = cms.double(3.0),
     GenMaxAbsEta     = cms.double(2.55),
     GenMaxVertR      = cms.double(3.5), # Maximum distance of particle production vertex from centre of CMS.
     GenMaxVertZ      = cms.double(25.0),
     GenPdgIds        = cms.vuint32(), # Only particles with these PDG codes used for efficiency measurement.

     # Additional cut on MC truth tracks used for algorithmic tracking efficiency measurements.
     # You should usually set this equal to value of L1TrackDef.MinStubLayers below, unless L1TrackDef.MinPtToReduceLayers
     # is < 10000, in which case, set it equal to (L1TrackDef.MinStubLayers - 1).
     GenMinStubLayers = cms.uint32(5)
  ),

  #=== Cuts applied to stubs before arriving in L1 track finding board.

  StubCuts = cms.PSet(
     # Reduce number of bits used by front-end chips to store stub bend info. (CMS is considering this option proposed by Seb. Viret).
     BendResReduced = cms.bool(True),
     # Don't use stubs with eta beyond this cut, since the tracker geometry makes it impossible to reconstruct tracks with them.
     MaxStubEta     = cms.double(2.4),
     # Don't use stubs whose measured Pt from bend info is significantly below HTArraySpec.HoughMinPt, where "significantly" means allowing for resolution in q/Pt derived from stub bend resolution specified below.
     KillLowPtStubs = cms.bool(True),
     # Bend resolution assumed by bend filter in units of strip pitch. Also used when assigning stubs to sectors if EtaPhiSectors.CalcPhiTrkRes=True. And by the bend filter if HTFillingRphi.UseBendFilter=True.
     # Suggested value: 1.19 if BendResReduced = false, or 1.30 if it is true.
     BendResolution = cms.double(1.25),
     # Additional contribution to bend resolution from its encoding into a reduced number of bits.
     # This number is the assumed resolution relative to the naive guess of its value.
     # It is ignored in BendResReduced = False.
     BendResolutionExtra = cms.double(0.0)
  ),

  #=== Optional Stub digitization.

  StubDigitize = cms.PSet(
     EnableDigitize  = cms.bool(False),  # Digitize stub coords? If not, use floating point coords.
     FirmwareType    = cms.uint32(1),    # 0 = Old Thomas 2-cbin data format, 1 = new Thomas data format used for daisy chain, 2-4 = reserved for demonstrator use, 9 = Systolic array data format.
     #
     #--- Parameters available in MP board.
     #
     PhiSectorBits   = cms.uint32(6),    # Bits used to store phi sector number
     PhiSBits        = cms.uint32(13),   # Bits used to store phiS coord.
     PhiSRange       = cms.double(0.3926990817),   # Range phiS coord. covers in radians.
     RtBits          = cms.uint32(10),   # Bits used to store Rt coord.
     RtRange         = cms.double(103.0382), # Range Rt coord. covers in units of cm.
     ZBits           = cms.uint32(12),   # Bits used to store z coord.
     ZRange          = cms.double(640.), # Range z coord. covers in units of cm.
     # The following four parameters do not need to be specified if FirmwareType = 1 (i.e., daisy-chain firmware) 
     DPhiBits        = cms.untracked.uint32(8),    # Bits used to store Delta(phi) track angle.
     DPhiRange       = cms.untracked.double(1.),   # Range Delta(phi) covers in radians.
     RhoBits         = cms.untracked.uint32(6),    # Bits used to store rho parameter.
     RhoRange        = cms.untracked.double(0.25), # Range rho parameter covers.
     #
     #--- Parameters available in GP board (excluding any in common with MP specified above).
     #
     PhiOBits        = cms.uint32(15),      # Bits used to store PhiO parameter.
     PhiORange       = cms.double(1.5707963268), # Range PhiO parameter covers.
     BendBits        = cms.uint32(6)        # Bits used to store stub bend.
  ),

  #=== Division of Tracker into phi sectors.

  PhiSectors = cms.PSet(
     NumPhiSectors      = cms.uint32(32),  # IRT - 32 or 64 are reasonable choices.
     ChosenRofPhi       = cms.double(58.), # Use phi of track at this radius for assignment of stubs to phi sectors & also for one of the axes of the r-phi HT. If ChosenRofPhi=0, then use track phi0.
     #--- You can set one or both the following parameters to True.
     UseStubPhi         = cms.bool(True),  # Require stub phi to be consistent with track of Pt > HTArraySpec.HoughMinPt that crosses HT phi axis?
     UseStubPhiTrk      = cms.bool(True),  # Require stub phi0 (or phi65 etc.) as estimated from stub bend, to lie within HT phi axis, allowing tolerance(s) specified below?
     AssumedPhiTrkRes   = cms.double(0.5), # Tolerance in stub phi0 (or phi65) assumed to be this fraction of phi sector width. (N.B. If > 0.5, then stubs can be shared by more than 2 phi sectors).
     CalcPhiTrkRes      = cms.bool(True),  # If true, tolerance in stub phi0 (or phi65 etc.) will be reduced below AssumedPhiTrkRes if stub bend resolution specified in StubCuts.BendResolution suggests it is safe to do so.
     HandleStripsPhiSec = cms.bool(False)  # If True, adjust algorithm to allow for uncertainty in stub (r,z) coordinate caused by length of 2S module strips when assigning stubs to phi sectors.
  ),

  #=== Division of Tracker into eta sectors

  EtaSectors = cms.PSet(
     # Optimized choice for 9 sectors - you MUST set chosenRofZ = 45 to use this.
     EtaRegions = cms.vdouble(-2.4, -2.0, -1.53, -0.98, -0.37, 0.37, 0.98, 1.53, 2.0, 2.4),
     # Mark's original choice of 5 sectors - use chosenRofZ = 65 with this
     #EtaRegions  = cms.vdouble(-2.4, -1.45, -0.61, 0.61, 1.45, 2.4),
     # Optimized choice for 5 sectors - use chosenRofZ = 65 with this
     #EtaRegions  = cms.vdouble(-2.4, -1.74, -0.71, 0.71, 1.74, 2.4),
     # Optimized choice for 7 sectors - use chosenRofZ = 65 with this
     #EtaRegions  = cms.vdouble(-2.4, -1.9, -1.3, -0.5, 0.5, 1.3, 1.9, 2.4),
     ChosenRofZ  = cms.double(45.),       # Use z of track at this radius for assignment of tracks to eta sectors & also for one of the axes of the r-z HT. Do not set to zero!
     BeamWindowZ = cms.double(15),        # Half-width of window assumed to contain beam-spot in z.
     HandleStripsEtaSec = cms.bool(False) # If True, adjust algorithm to allow for uncertainty in stub (r,z) coordinate caused by length of 2S module strips when assigning stubs to eta sectors.
  ),

  #=== r-phi Hough transform array specifications.

  HTArraySpecRphi = cms.PSet(
     HoughMinPt      = cms.double(3.0), # Min track Pt that Hough Transform must find. Also used by StubCuts.KillLowPtStubs and by EtaPhiSectors.UseStubPhi.
     HoughNbinsPt    = cms.uint32(32),  # HT array dimension in track q/Pt. Ignored if HoughNcellsRphi > 0.
     HoughNbinsPhi   = cms.uint32(32),  # HT array dimension in track phi0 (or phi65 or any other track phi angle. Ignored if HoughNcellsRphi > 0.
     HoughNcellsRphi = cms.int32(-1),   # If > 0, then parameters HoughNbinsPt and HoughNbinsPhi will be calculated from the constraints that their product should equal HoughNcellsRphi and their ratio should make the maximum |gradient|" of stub lines in the HT array equal to 1. If <= 0, then HoughNbinsPt and HoughNbinsPhi will be taken from the values configured above.
     EnableMerge2x2  = cms.bool(False), # Groups of neighbouring 2x2 cells in HT will be treated as if they are a single large cell? N.B. You can only enable this option if your HT array has even numbers of bins in both dimensions. 
     MaxPtToMerge2x2 = cms.double(6.),  # but only cells with pt < MaxPtToMerge2x2 will be merged in this way (irrelevant if EnableMerge2x2 = false).
     NumSubSecsEta   = cms.uint32(1)    # Subdivide each sector into this number of subsectors in eta within r-phi HT.
  ),

  #=== r-z Hough transform array specifications.

  HTArraySpecRz = cms.PSet(
     EnableRzHT    = cms.bool(False),  # If true, find tracks with both r-phi & r-z HTs. If false, use only r-phi HT. If false, other parameters in this sector irrelevant.
     HoughNbinsZ0  = cms.uint32(0),   # HT array dimension in track z0. Ignored if HoughNcellsRz > 0.
     HoughNbinsZ65 = cms.uint32(0),   # HT array dimension in track z65 (or any other r-z track related variable). Ignored if HoughNcellsRz > 0.
     HoughNcellsRz = cms.int32(1024)  # If > 0, then parameters HoughNbinsZ0 and HoughNbinsZ65 will be calculated from the constraints that their product should equal HoughNcellsRz and their ratio should make the maximum |gradient|" of stub lines in the HT array equal to 1. If <= 0, then HoughNbinsZ0 and HoughNbinsRz will be taken from the values configured above.
  ),

  #=== Rules governing how stubs are filled into the r-phi Hough Transform array.

  HTFillingRphi = cms.PSet(
     # If True, adjust algorithm to allow for uncertainty in stub (r,z) coordinate caused by length of 2S module strips when filling r-phi HT with stubs.
     HandleStripsRphiHT   = cms.bool(False),
     # Take all cells in r-phi HT array crossed by line corresponding to each stub (= 0) or take only some to reduce rate at cost
     # of efficiency ( > 0). If this option is > 0, it can be 1 or 2, corresponding to different algorithms for rejecting
     # some of the cells. "1" is an algorithm invented by Ian, whereas "2" corresponds to Thomas' 1st firmware implementation which only handled 1 cell per HT column.
     # Suggest setting KillSomeHTCellsRphi=1 (=0) if HTArraySpec.ChosenRofPhi=0 (>0)
     KillSomeHTCellsRphi  = cms.uint32(0),
     # Use filter in each r-phi HT cell, filling it only with stubs that have consistent bend information?
     # The assumed bend resolution is specified in StubCuts.BendResolution.
     UseBendFilter        = cms.bool(True),
     # A filter is used each HT cell, which prevents more than the specified number of stubs being stored in the cell. (Reflecting memory limit of hardware).
     MaxStubsInCell       = cms.uint32(99999), # Setting this to anything more than 99 disables this option
     #MaxStubsInCell      = cms.uint32(16),    # set it equal to value used in hardware.
     # If BusySectorKill = True, and more than BusySectorNumStubs stubs are assigned to tracks by an r-phi HT array, then the excess tracks are killed, with lowest Pt ones killed first. This is because hardware has finite readout time.
     BusySectorKill       = cms.bool(False),
     BusySectorNumStubs   = cms.uint32(210),
     # If this is True, then the BusySectorNumStubs cut is applied to +ve and -ve charge track seperately. (Irrelevant if BusySectorKill = False).
     BusySectorEachCharge = cms.bool(False)
  ),

  #=== Rules governing how stubs are filled into the r-z Hough Transform array. (Irrelevant if HTArraySpecRz.enableRzHT = false)

  HTFillingRz = cms.PSet(
     # If True, adjust algorithm to allow for uncertainty in stub (r,z) coordinate caused by length of 2S module strips when filling r-z HT with stubs.
     HandleStripsRzHT  = cms.bool(True),
     # Take all cells in r-z HT array crossed by line corresponding to each stub (= 0) or take only some to reduce rate at cost
     # of efficiency ( > 0). If this option is > 0, it can be 1 or 2, corresponding to different algorithms for rejecting
     # some of the cells.
     KillSomeHTCellsRz = cms.uint32(0)
  ),

  #=== Options controlling r-z track filters (or any other track filters run after the Hough transform, as opposed to inside it).

  RZfilterOpts = cms.PSet(
     # Use filter in each r-phi HT cell, filling it only with stubs that have consistent rapidity?
     UseEtaFilter        = cms.bool(False),
     # Use filter in each HT cell using only stubs which have consistent value of ZTrk
     UseZTrkFilter       = cms.bool(False),
     # Filter Stubs in each HT cell using a tracklet seed algorithm
     UseSeedFilter       = cms.bool(True),
     #--- Options for Ztrk filter, (so only relevant if UseZtrkFilter=true).
     # Use z of track at this radius for ZTrkFilter. 
     ChosenRofZFilter    = cms.double(23.),
     #--- Options relevant for Seed filter, (so only relevant if useSeedFilter=true).
     # Added resolution for a tracklet-like filter algorithm, beyond that estimated from hit resolution. 
     SeedResolution      = cms.double(0.),
     # Store stubs compatible with all possible good seed.
     KeepAllSeed         = cms.bool(True),
     # Maximum number of seed combinations to bother checking per track candidate.
     MaxSeedCombinations = cms.uint32(13),
     # Reject tracks whose estimated rapidity from seed filter is inconsistent range of with eta sector. (Kills some duplicate tracks).
     zTrkSectorCheck     = cms.bool(True)
  ),

  #=== Rules for deciding when the track finding has found an L1 track candidate

  L1TrackDef = cms.PSet(
     # Min. number of layers in HT cell that must have stubs to produce a reco track candidate.
     MinStubLayers        = cms.uint32(5),
     # Reduce MinStubLayers by 1 for HT cells corresponding to track pt above this cut.
     # If this is set to >= 10000., this option is disabled.
     MinPtToReduceLayers  = cms.double(99999.),
     # Define layers using layer ID (true) or by bins in radius of 5 cm width (false).
     UseLayerID           = cms.bool(True),
     # Reduce this layer ID, so that it takes no more than 8 different values in any eta region (simplifies firmware).
     ReducedLayerID       = cms.bool(True)
  ),

  #=== Specification of algorithm to eliminate duplicate tracks.

  DupTrkRemoval = cms.PSet(
    #--- Specify which duplicate removal algorithm(s) to run:  option 0 means disable duplicate track removal, whilst > 0 runs a specific algorithm.
    # Algorithm used for duplicate removal of 2D tracks produced by r-phi HT. Assumed to run before tracks are output from HT board.
    DupTrkAlgRphi = cms.uint32(0),
    #DupTrkAlgRphi = cms.uint32(12),
    # Algorithm used for duplicate removal run on 2D tracks produced by r-z HT corresponding to a single r-phi HT track. Assumed to run before tracks are output from HT board.
    DupTrkAlgRz   = cms.uint32(0),
    # Algorithm run on all 3D tracks within each sector after r-z HT or r-z seed filter.
    #DupTrkAlgRzSeg = cms.uint32(0),
    DupTrkAlgRzSeg = cms.uint32(8),
    # Algorithm run on tracks after the track helix fit has been done.
    DupTrkAlgFit   = cms.uint32(0),
    #DupTrkAlgFit   = cms.uint32(50),
    #--- Options used by individual algorithms.
    # Parameter for OSU duplicate-removal algorithm
    # Specifies minimum number of independent stubs to keep candidate in comparison in Algo 3
    DupTrkMinIndependent = cms.uint32(3),
    # Parameter for "inverse" OSU duplicate-removal algorithm
    # Specifies minimum number of common stubs in same number of layers to keep smaller candidate in comparison in Algos 5-9 + 15-16.
    DupTrkMinCommonHitsLayers = cms.uint32(5),
    # Reduced ChiSq cut in linear RZ fit in Algo 4
    DupTrkChiSqCut = cms.double(99999.),
    # Max diff in qOverPt of 2 tracks in Algo15
    DupMaxQOverPtScan = cms.double(0.025),
    # Max diff in phi0 range of 2 tracks in Algo15
    DupMaxPhi0Scan = cms.double(0.01),
    # Max diff in z0 of 2 tracks in Algo15
    DupMaxZ0Scan = cms.double(0.2),
    # Max diff in tanLambda of 2 tracks in Algo 15
    DupMaxTanLambdaScan = cms.double(0.01)
  ),

  #=== Rules for deciding when a reconstructed L1 track matches a MC truth particle (i.e. tracking particle).

  TrackMatchDef = cms.PSet(
     #--- Three different ways to define if a tracking particle matches a reco track candidate. (Usually, set two of them to ultra loose).
     # Min. fraction of matched stubs relative to number of stubs on reco track.
     MinFracMatchStubsOnReco  = cms.double(-99.),
     # Min. fraction of matched stubs relative to number of stubs on tracking particle.
     MinFracMatchStubsOnTP    = cms.double(-99.),
     # Min. number of matched layers.
     MinNumMatchLayers        = cms.uint32(5),
     # Associate stub to TP only if the TP contributed to both its clusters? (If False, then associate even if only one cluster was made by TP).
     StubMatchStrict          = cms.bool(False)
  ),

  #=== Track Fitting Algorithm Settings.

  TrackFitSettings = cms.PSet(
     #
     #--- Options applicable to all track fitters ---
     #
     # Track Fitting algortihms to use. You can run several in parallel.
     # (TrackFitLinearAlgo & ChiSquared* are chi2 fits, KF* is a Kalman filter fit, and LinearRegression is simplified 4 
     # parameter fit that neglects the hit uncertainties. The number 4 or 5 indicates if 4 or 5 helix parameters are fitted).
     # WARNING: KF5ParamsComb crashes, so don't use it!
     TrackFitters = cms.vstring(
                                "TrackFitLinearAlgo4",
                                #"TrackFitLinearAlgo5",
                                #"ChiSquared4ParamsApprox",
                                #"ChiSquared5ParamsApprox",
                                #"ChiSquared4ParamsTrackletStyle",
                                "KF4ParamsComb",
                                #"KF5ParamsComb",
                                "LinearRegression"
                                ),
     # Cut on chi2/dof of fitted track when making histograms.
     Chi2OverNdfCut = cms.double(999999.),
     # Print detailed summary of track fit performance at end of job (as opposed to a brief one). 
     DetailedFitOutput = cms.bool(False),
     #
     #--- Options for chi2 & Linear Regression track fitters ---
     #
     # Number of fitting iterations to undertake. (15 is not realistic in hardware, but is necessary to kill bad hits)
     NumTrackFitIterations = cms.uint32(15),
     # Optionally kill hit with biggest residuals in track fit (occurs after the first fit, so three iterations would have two killings). (Only used by chi2 track fit).
     KillTrackFitWorstHit  = cms.bool(True),
     # Cuts in standard deviations used to kill hits with big residuals during fit. If the residual exceeds the "General" cut, the hit is killed providing it leaves the track with enough hits to survive. If the residual exceeds the "Killing" cut, the hit is killed even if that kills the track.
     GeneralResidualCut = cms.double(3.0),
     KillingResidualCut = cms.double(20.0),
     #
     #--- Additional options for Linear Regression track fitter ---
     #
     # If False: residual of a stub is the max of its r-phi & r-z residuals. 
     # If True: the residual is the mean of these residuals.
     CombineResiduals                = cms.bool( True ),
     # Correct stub phi coordinate for higher orders in circle expansion, so that a trajectory is straight in r-phi.
     LineariseStubPosition           = cms.bool( True ),
     # Checks if the fitted track is consistent with the sector, if not it will be not accepted.
     CheckSectorConsistency          = cms.bool( False ),
     # Checks if the fitted track r phi parameter  are consistent with the HT candidate parameter within in range of +- 2 cells.
     CheckHTCellConsistency          = cms.bool( False ),
     # Tracks must have stubs in at least this number of PS layers.
     MinPSLayers                     = cms.uint32( 2 ),
     FitPerfectCandidatesOnly        = cms.bool( False ),
     #
     #--- Options for Kalman filter track fitters ---
     #
     # larger number has more debugging outputs.
     KalmanDebugLevel                = cms.uint32(0),
     # Internal histograms are filled if it is True
     KalmanFillInternalHists         = cms.bool(True),
     # Multiple scattering factor.  Not working. Set to 0.
     KalmanMultipleScatteringFactor  = cms.double(0.0),
     # A stub which is inconsistent with the state is not processed for the state. Cut value on chisquare from the forcast is set.
     KalmanValidationGateCutValue    = cms.double(150),
     # Best candidate is selected from the candidates with the more number of stubs if this is Ture. Chi2 only selection is better to remove combinatorial hits.
     KalmanSelectMostNumStubState    = cms.bool(False),
     # After this number of stubs in the state, the rest of the stubs are not processed and go to the next state.
     KalmanMaxNumNextStubs           = cms.uint32(999),
     # Allowed number of virtual stubs. (CAUTION: kalmanState::nVirtualStubs() return the # of virtual stubs including the seed. The number on this option does not include the seed.)
     KalmanMaxNumVirtualStubs        = cms.uint32(1),
     # The number of states for each virtual stub list is restricted to this number.  The lists are truncated to this number of states.
     KalmanMaxNumStatesCutValue      = cms.uint32(10),
     # The state is removed if the reduced chisquare is more than this number.
     KalmanStateReducedChi2CutValue  = cms.double(100)
  ),

  # Debug printout
  Debug  = cms.uint32(0) #(0=none, 1=print tracks/sec, 2=show filled cells in HT array in each sector of each event, 3=print all HT cells each TP is found in, to look for duplicates, 4=print missed tracking particles by r-z filters, 5 = show debug info about duplicate track removal, 6 = show debug info about fitters)
)
