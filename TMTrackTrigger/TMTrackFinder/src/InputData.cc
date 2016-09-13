#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"

#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"

// TTStubAssociationMap.h forgets to two needed files, so must include them here ...
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"

#include "TMTrackTrigger/TMTrackFinder/interface/InputData.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"

#include <map>

using namespace std;
 
InputData::InputData(const edm::Event& iEvent, const edm::EventSetup& iSetup, Settings* settings) {

  vTPs_.reserve(2500);
  vStubs_.reserve(35000);
  vAllStubs_.reserve(35000);

  // Get TrackingParticle info

  edm::Handle<TrackingParticleCollection> tpHandle;
  iEvent.getByLabel("mix", "MergedTrackTruth", tpHandle );

  unsigned int tpCount = 0;
  for (unsigned int i = 0; i < tpHandle->size(); i++) {
    TrackingParticlePtr tpPtr(tpHandle, i);
    // Store the TrackingParticle info, using class TP to provide easy access to the most useful info.
    TP tp(tpPtr, tpCount, settings);
    // Only bother storing tp if it could be useful for tracking efficiency or fake rate measurements.
    if (tp.use()) {
      vTPs_.push_back( tp );
      tpCount++;
    }
  }

  // Also create map relating edm::Ptr<TrackingParticle> to TP.

  map<edm::Ptr< TrackingParticle >, const TP* > translateTP;

  for (const TP& tp : vTPs_) {
    TrackingParticlePtr tpPtr(tp);
    translateTP[tpPtr] = &tp;
  }

  // Get the tracker geometry info needed to unpack the stub info.

  edm::ESHandle<StackedTrackerGeometry> stackedGeometryHandle;
  iSetup.get<StackedTrackerGeometryRecord>().get( stackedGeometryHandle );

  const StackedTrackerGeometry*  stackedGeometry = stackedGeometryHandle.product();

  // Get stub info, by looping over modules and then stubs inside each module.
  // Also get the association map from stubs to tracking particles.

  edm::Handle<DetSetVec>       ttStubHandle;
  edm::Handle<TTStubAssMap>    mcTruthTTStubHandle;
  edm::Handle<TTClusterAssMap> mcTruthTTClusterHandle;
  iEvent.getByLabel("TTStubsFromPixelDigis"            , "StubAccepted"    , ttStubHandle           );
  iEvent.getByLabel("TTStubAssociatorFromPixelDigis"   , "StubAccepted"    , mcTruthTTStubHandle    );
  iEvent.getByLabel("TTClusterAssociatorFromPixelDigis", "ClusterAccepted" , mcTruthTTClusterHandle );

  unsigned int stubCount = 0;
  for (DetSetVec::const_iterator p_module = ttStubHandle->begin(); p_module != ttStubHandle->end(); p_module++) {
    for (DetSet::const_iterator p_ttstub = p_module->begin(); p_ttstub != p_module->end(); p_ttstub++) {
      TTStubRef ttStubRef = edmNew::makeRefTo(ttStubHandle, p_ttstub );
      // Store the Stub info, using class Stub to provide easy access to the most useful info.
      Stub stub(ttStubRef, stubCount, settings, stackedGeometry);
      // Also fill truth associating stubs to tracking particles.
      //      stub.fillTruth(vTPs_, mcTruthTTStubHandle, mcTruthTTClusterHandle); 
      stub.fillTruth(translateTP, mcTruthTTStubHandle, mcTruthTTClusterHandle); 
      vAllStubs_.push_back( stub );
      stubCount++;
    }
  }

  // Produced reduced list containing only the subset of stubs that the user has declared will be 
  // output by the front-end readout electronics.
  for (const Stub& s : vAllStubs_) {
    if (s.frontendPass()) vStubs_.push_back( &s );
  }


  // Note list of stubs produced by each tracking particle.

  // (By passing vAllStubs_ here instead of vStubs_, it means that any algorithmic efficiencies
  // measured will be reduced if the tightened frontend electronics cuts, specified in section StubCuts
  // of Analyze_Defaults_cfi.py, are not 100% efficient).
  for (unsigned int j = 0; j < vTPs_.size(); j++) {
    vTPs_[j].fillTruth(vAllStubs_);
  }
}
