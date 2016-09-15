#ifndef __TMTRACKPRODUCER_H__
#define __TMTRACKPRODUCER_H__

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/Demonstrator/interface/HardwareStub.h"
#include "DataFormats/Demonstrator/interface/HardwareTrack.h"

#include <vector>
#include <map>
#include <string>


class Settings;
class Histos;
class TrackFitGeneric;

class TMTrackProducer : public edm::EDProducer {

public:
  explicit TMTrackProducer(const edm::ParameterSet&);	
  ~TMTrackProducer(){}

private:

/*CMSSW_8_MIGRATION*/ //   typedef std::vector< TTTrack< Ref_PixelDigi_ > > TTTrackCollection;
  typedef std::vector<l1t::HardwareStub>           HwStubCollection;
  typedef std::vector<l1t::HardwareTrack>          HwTrackCollection;

  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

private:

  Settings *settings_;
  Histos   *hists_;
  std::map<std::string, TrackFitGeneric*> fitterWorkerMap_;
};
#endif

