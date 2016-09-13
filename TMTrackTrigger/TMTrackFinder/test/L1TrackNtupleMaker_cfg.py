#########################################################################################################
# This is a copy of Louise Skinnari's script L1TrackNupleMaker_cfg.py obtained via
# getLouiseAnalysisCode.csh .
# It has been modified to read our MC and to run our L1 track producer.
#
# To run execute do
# cmsRun L1TrackNtupleMaker_cfg.py Events=50 inputMC=Samples/Muons/PU0.txt histFile=outputHistFile.root outEdmFile=outputEdm.root trkFitAlgo=TrackFitLinearAlgo4
# where the arguments take default values if you don't specify them. You can change defaults below.
#
# The option trkFitAlgo should be equal to the name of the TMT track fitting algorithm that you wish to use,
# which must be one of the fitters specified by option TrackFitters in TMTrackProducer_Defaults_cfi.py.
# Alternatively, you can set it equal to "Tracklet4" or "Tracklet5" to use the Tracklet groups 
# tracks with 4 or 5 helix parameters.
#########################################################################################################


############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
import os

process = cms.Process("L1TrackNtuple")
 
#################################################################################
# Options such as choice of MC sample (specific to us)
################################################################################

options = VarParsing.VarParsing ('analysis')

#--- Specify input MC
#options.register('inputMC', '../../../SamplesCMSseb/Muon_Pt100/PU140.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
#options.register('inputMC', '../../../SamplesCMSseb/Electron_Pt35/PU140.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
#options.register('inputMC', '../../../SamplesCMSseb/StubFix/TTbar/PU140.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
options.register('inputMC', '../../../SamplesCMS/StubFix/TTbar/PU140.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")

#--- Specify number of events to process.
options.register('Events',100,VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,"Number of Events to analyze")

#--- Specify name of output histogram file.
options.register('histFile','Hist.root',VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"Name of output histogram file")

#--- Specify name of output EDM file contains stubs associated with tracks to compare with hardware.
#--- If the name is equal to a null string, no EDM file will be written.
options.register('outEdmFile','',VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"Name of output EDM file")
#options.register('outEdmFile','outputEdm.root',VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"Name of output EDM file")

options.register('trkFitAlgo','TrackFitLinearAlgo4',VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"Nuame of track helix fit algorithm")

options.parseArguments()

#--- Decode which track fitting algorithm results are wanted for

# Note name of track fitting algorithm.
trkFitAlgo = options.trkFitAlgo
# Are results wanted for the 4 or 5 parameter helix fit?
if ( (trkFitAlgo.find("4") == -1 and trkFitAlgo.find("5") != -1) or (trkFitAlgo.find("5") == -1 and trkFitAlgo.find("4") != -1) ) :
    # If this is true, it was a 5 param fit; else 4 param.
    fit5param = (trkFitAlgo.find("5") != -1)
else:
    print "### ERROR: Unable to figure out if helix fit used 4 or 5 params from its algorithm name"
# Are our TMT tracks to be used? If not, then the Tracklet group's tracks will be used.
useTMTtracks = (trkFitAlgo.find("Tracklet") == -1)
# Figure out name of our TTTrack collection
if (useTMTtracks):
  trackName = "TML1Tracks" + trkFitAlgo
else:
  trackName = "Level1TTTracks"

if (fit5param):
    print "### Will produce results for TTTrack collection ",trackName," with 5 helix parameters."
else:
    print "### Will produce results for TTTrack collection ",trackName," with 4 helix parameters."

############################################################
# input and output (specific to us, as we use our own MC)
############################################################

list = FileUtils.loadListFromFile(options.inputMC)
readFiles = cms.untracked.vstring(*list)
# readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.Events) )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.histFile)
)

process.source = cms.Source ("PoolSource",
                            fileNames = readFiles,
                            secondaryFileNames = secFiles
                            )

process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))
 
############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

############################################################
# Path definitions & schedule
############################################################

#run the tracking
BeamSpotFromSim = cms.EDProducer("BeamSpotFromSimProducer")
process.TT_step = cms.Path(process.TrackTriggerTTTracks)
process.TTAssociator_step = cms.Path(process.TrackTriggerAssociatorTracks)

############################################################
# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# e.g. single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13 
#      pions in jets = 6
#      taus = 15
#      all TPs = 1
############################################################

# Modified this to specify our tracks & associator to be used.
process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(False),      # printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   # save all L1 tracks, not just truth matched to primary particle
                                       SaveStubs = cms.bool(False),      # save some info for *all* stubs
                                       L1Tk_nPar = cms.int32(4),         # use 4 or 5-parameter L1 track fit ??
                                       L1Tk_minNStub = cms.int32(5),     # L1 tracks with >= 4 stubs
                                       TP_minNStub = cms.int32(4),       # require TP to have >= minNstub (efficiency denominator)
                                       TP_minNStubLayer = cms.int32(4),  # require TP to have stubs in >= X layers/disks
                                       TP_minPt = cms.double(1.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.4),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       L1TrackInputTag = cms.InputTag("TTTracksFromPixelDigis", "Level1TTTracks" ),               # TTTrack input
                                       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), # MCTruth input 
                                       )

process.ana = cms.Path(process.L1TrackNtuple)

# Specific to us:
# Optional output of EDM file containing stubs associated to tracks to compare with hardware.
# Can use this to store edm-file instead of root-ntuple.
if (options.outEdmFile != "") :

    process.out = cms.OutputModule( "PoolOutputModule",
                                    fileName = cms.untracked.string(options.outEdmFile),
                                    fastCloning = cms.untracked.bool( False ),
                                    outputCommands = cms.untracked.vstring('drop *',
                                                                           'keep *_*_Level1PixelTracks_*',
                                                                           'keep *_*_Level1TTTracks_*',
                                                                           'keep *_*_StubAccepted_*',
                                                                           'keep *_*_ClusterAccepted_*',
                                                                           'keep *_*_MergedTrackTruth_*',
                                                                           'keep *_TMTrackProducer_*_*')
    )

    process.FEVToutput_step = cms.EndPath(process.out)

############################################################
# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
############################################################

# OLD, FOR TTI SAMPLES
# from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023TTI
# process = cust_2023TTI(process)

# # NEW, FOR SCOPE DOCUMENT REL VAL SAMPLES
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D
process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)

#=== EVERYTHING BELOW WAS ADDED BY THE TMT GROUP.

# Set number of helix fit parameters that analysis code should get results for.
if (fit5param):
    process.L1TrackNtuple.L1Tk_nPar = cms.int32(5)


# If use wants plots for our TMT L1 tracks instead of Tracklet group's ones, then the following is run.
if (useTMTtracks):
    #--- Load code that produces our L1 tracks and makes corresponding histograms.
    #--- Either use this one for studies of the final 2025 system.
    process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_cff')
    #--- Or this one for studies of the 2015 demonstrator with Thomas Schuh's "daisy chain" firmware.
    #process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_Demo_DaisyChain_cff')
    #--- Or this one for studies of the 2015 demonstrator with Thomas Schuh's "2-c-bin" firmware.
    #process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_Demo_Thomas_cff')
    #--- Or this one for studies of the 2015 demonstrator with systolic array firmware.
    #process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_Demo_Systolic_cff')

    #--- Add boolean, indicating if output EDM file will be written, to cfg params, so is available to C++.
    process.TMTrackProducer.WriteOutEdmFile = cms.untracked.bool( (options.outEdmFile != "") )

    #--- Optionally override default configuration parameters here (example given of how).

    #process.TMTrackProducer.HTArraySpecRz.EnableRzHT = cms.bool(True)

    # Change tracking step to run our L1 track reconstruction algorithm instead of the tracklet group's one.
    process.TT_step = cms.Path(process.BeamSpotFromSim+process.TMTrackProducer)

    # Change track association to truth step to use our tracks instead of tracklet group's ones.
    process.TMTrackAssociator = process.TTTrackAssociatorFromPixelDigis.clone(
        TTTracks = cms.VInputTag(cms.InputTag("TMTrackProducer", trackName))
    )
    process.TTAssociator_step = cms.Path(process.TMTrackAssociator)

    # Run analysis code over our L1 tracks.
    process.L1TrackNtuple.L1TrackInputTag      = cms.InputTag("TMTrackProducer", trackName)
    # And use the associator that matches our tracks to the truth tracks.
    process.L1TrackNtuple.MCTruthTrackInputTag = cms.InputTag("TMTrackAssociator", trackName)

process.schedule = cms.Schedule(process.TT_step,process.TTAssociator_step,process.ana)
