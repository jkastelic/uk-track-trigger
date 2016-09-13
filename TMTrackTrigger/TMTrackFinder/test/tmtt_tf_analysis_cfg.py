#################################################################################################
# To run execute do
# cmsRun tmtt_tf_analysis_cfg.py Events=50 inputMC=Samples/Muons/PU0.txt histFile=outputHistFile.root outEdmFile=outputEdm.root
# where the arguments take default values if you don't specify them. You can change defaults below.
#################################################################################################

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

process.load('Configuration.Geometry.GeometryExtended2023Pixel_cff')
process.load('Configuration.Geometry.GeometryExtended2023PixelReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')

process.load("FWCore.MessageLogger.MessageLogger_cfi")

options = VarParsing.VarParsing ('analysis')

#--- Specify input MC
#options.register('inputMC', '../../../Samples/Muons_fixed_pT_10GeV/PU0.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
#options.register('inputMC', '../../../Samples/Electrons_FixedPt_10GeV/PU0.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
#options.register('inputMC', '../../../Samples/TTbar/StubFix/PU140.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")
options.register('inputMC', '../../../SamplesCMS/StubFix/TTbar/PU140.txt', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, "Files to be processed")

#--- Specify number of events to process.
options.register('Events',100,VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,"Number of Events to analyze")

#--- Specify name of output histogram file.
options.register('histFile','Hist.root',VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"Name of output histogram file")

#--- Specify name of output EDM file contains stubs associated with tracks to compare with hardware.
#--- If the name is equal to a null string, no EDM file will be written.
options.register('outEdmFile','',VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"Name of output EDM file")
#options.register('outEdmFile','outputEdm.root',VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,"Name of output EDM file")

options.parseArguments()

#--- input and output

list = FileUtils.loadListFromFile(options.inputMC)
readFiles = cms.untracked.vstring(*list)
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

#--- Load code that produces our L1 tracks and makes corresponding histograms.

#--- Either use this one for studies of the final 2025 system.
process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_cff')
#--- Or this one for studies of the 2015 demonstrator with Thomas Schuh's "daisy chain" firmware.
#process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_Demo_DaisyChain_cff')
#--- Or this one for studies of the 2015 demonstrator with Thomas Schuh's old firmware.
#process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_Demo_Thomas_cff')
#--- Or this one for studies of the 2015 demonstrator with systolic array firmware.
#process.load('TMTrackTrigger.TMTrackFinder.TMTrackProducer_Demo_Systolic_cff')

#--- Add boolean, indicating if output EDM file will be written, to cfg params, so is available to C++.
process.TMTrackProducer.WriteOutEdmFile = cms.untracked.bool( (options.outEdmFile != "") )

#--- Optionally override default configuration parameters here (example given of how).

#process.TMTrackProducer.HTArraySpecRz.EnableRzHT = cms.bool(True)

process.p = cms.Path(process.TMTrackProducer)

# Optional output of EDM file containing stubs associated to tracks to compare with hardware.
if (options.outEdmFile != "") :

  process.out = cms.OutputModule("PoolOutputModule",
              fileName = cms.untracked.string(options.outEdmFile),
              outputCommands = cms.untracked.vstring('drop *',
                                                     'keep *_TMTrackProducer_*_*')
  )

  process.e = cms.EndPath(process.out)
