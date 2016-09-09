import FWCore.ParameterSet.Config as cms

#=== Import default values for all parameters.

from TMTrackTrigger.TMTrackFinder.TMTrackProducer_Defaults_cfi import TMTrackProducer

#=== Change below any parameters you don't want to take their default values.

#--- Use Thomas's improved Hough transform (high granularity in cells with Pt > 6 GeV; and 2 eta subsectors;
#--- and outputting +ve and -ve charged tracks on separate links.

#TMTrackProducer.HTArraySpecRphi.HoughNbinsPt = cms.uint32(64) 
#TMTrackProducer.HTArraySpecRphi.HoughNbinsPhi = cms.uint32(64)
#TMTrackProducer.HTArraySpecRphi.EnableMerge2x2  = cms.bool(True)
#TMTrackProducer.HTArraySpecRphi.MaxPtToMerge2x2 = cms.double(6.)

#TMTrackProducer.HTArraySpecRphi.NumSubSecsEta   = cms.uint32(2)

#TMTrackProducer.HTFillingRphi.BusySectorEachCharge = cms.bool(True)

#--- Kill tracks that HT doesn't have time to output during time-multiplexed period.
#--- (Reduces efficiency, but is more realistic).

#TMTrackProducer.HTFillingRphi.BusySectorKill       = cms.bool(True)

#--- Switch on Ian's duplicate track removal algorithm that runs after the track fitter.
#--- Requires all other duplicate track removal to be disabled.

#TMTrackProducer.DupTrkRemoval.DupTrkAlgRphi  = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgRz    = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgRzSeg = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgFit   = cms.uint32(50)

#--- Switch off parts of the track reconstruction chain.

#TMTrackProducer.RZfilterOpts.UseSeedFilter = cms.bool(False)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgRphi  = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgRz    = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgRzSeg = cms.uint32(0)
#TMTrackProducer.DupTrkRemoval.DupTrkAlgFit   = cms.uint32(0)
#TMTrackProducer.TrackFitSettings.TrackFitters = cms.vstring()



