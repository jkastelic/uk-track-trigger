#include <TMTrackTrigger/TMTrackFinder/interface/TMTrackProducer.h>
#include <TMTrackTrigger/TMTrackFinder/interface/InputData.h>
#include <TMTrackTrigger/TMTrackFinder/interface/Settings.h>
#include <TMTrackTrigger/TMTrackFinder/interface/Histos.h>
#include <TMTrackTrigger/TMTrackFinder/interface/Sector.h>
#include <TMTrackTrigger/TMTrackFinder/interface/HTpair.h>
#include <TMTrackTrigger/TMTrackFinder/interface/KillDupFitTrks.h>
#include <TMTrackTrigger/TMTrackFinder/interface/TrackFitGeneric.h>
#include <TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h>
#include <TMTrackTrigger/TMTrackFinder/interface/L1fittedTrk4and5.h>
#include <TMTrackTrigger/TMTrackFinder/interface/ConverterToTTTrack.h>
#include "TMTrackTrigger/TMTrackFinder/interface/HTcell.h"
#include "TMTrackTrigger/TMTrackFinder/interface/DemoOutput.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "boost/numeric/ublas/matrix.hpp"
#include <iostream>
#include <vector>
#include <set>

using namespace std;
using  boost::numeric::ublas::matrix;

TMTrackProducer::TMTrackProducer(const edm::ParameterSet& iConfig) {
  // Get configuration parameters
  settings_ = new Settings(iConfig);

  // Tame debug printout.
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(4);

  // Book histograms.
  hists_ = new Histos( settings_ );
  hists_->book();

  // Create track fitting algorithm (& internal histograms if it uses them)
  for (const string& fitterName : settings_->trackFitters()) {
    fitterWorkerMap_[ fitterName ] = TrackFitGeneric::create(fitterName, settings_);
    fitterWorkerMap_[ fitterName ]->bookHists(); 
  }

  //--- Define EDM output to be written to file (if required) 

  // L1 tracks found by Hough Transform without any track fit.
  produces< std::vector< TTTrack< Ref_PixelDigi_ > > >( "TML1TracksHT" ).setBranchAlias("TML1TracksHT");
  // L1 tracks after track fit by each of the fitting algorithms under study
  for (const string& fitterName : settings_->trackFitters()) {
    string edmName = string("TML1Tracks") + fitterName;
    produces< std::vector< TTTrack< Ref_PixelDigi_ > > >(edmName).setBranchAlias(edmName);
  }
  // Stubs for used by comparison software, to compare hardware with software.
  produces<HwStubCollection>("OutputSimStub").setBranchAlias("OutputSimStubs");
  produces<HwStubCollection>("AllOutputSimStub").setBranchAlias("AllOutputSimStubs");
  produces<HwTrackCollection>("EfficiencyTrack").setBranchAlias("EfficiencyTracks");
  produces<HwTrackCollection>("AlgoEfficiencyTrack").setBranchAlias("AlgoEfficiencyTracks");
}


void TMTrackProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) 
{
  // Get the B-field and store its value in the Settings class.

  edm::ESHandle<MagneticField> magneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
  const MagneticField* theMagneticField = magneticFieldHandle.product();
  float bField = theMagneticField->inTesla(GlobalPoint(0,0,0)).z(); // B field in Tesla.
  cout<<endl<<"--- B field = "<<bField<<" Tesla ---"<<endl<<endl;

  settings_->setBfield(bField);

  // Initialize track fitting algorithm at start of run (especially with B-field dependent variables).
  for (const string& fitterName : settings_->trackFitters()) {
    fitterWorkerMap_[ fitterName ]->initRun(); 
  }
}

void TMTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Note useful info about MC truth particles and about reconstructed stubs .
  InputData inputData(iEvent, iSetup, settings_);

  const vector<TP>&          vTPs   = inputData.getTPs();
  const vector<const Stub*>& vStubs = inputData.getStubs(); 

  cout<<"INPUT #TPs = "<<vTPs.size()<<" #STUBs = "<<vStubs.size()<<endl;

  //=== Fill histograms with stubs and tracking particles from input data.
  hists_->fillInputData(inputData);

  // Creates matrix of Sector objects, which decide which stubs are in which (eta,phi) sector
  matrix<Sector>  mSectors(settings_->numPhiSectors(), settings_->numEtaRegions());
  // Create matrix of Hough-Transform arrays, with one-to-one correspondence to sectors.
  matrix<HTpair>  mHtPairs(settings_->numPhiSectors(), settings_->numEtaRegions());

  //=== Initialization
  // Create utility for converting L1 tracks from our private format to official CMSSW EDM format.
  const ConverterToTTTrack converter(settings_);
  // Storage for EDM L1 track collection to be produced from Hough transform output (no fit).
  std::auto_ptr<TTTrackCollection>  htTTTracksForOutput(new TTTrackCollection);
  // Storage for EDM L1 track collection to be produced from fitted tracks (one for each fit algorithm being used).
  // auto_ptr cant be stored in std containers, so use C one, together with map noting which element corresponds to which algorithm.
  const unsigned int nFitAlgs = settings_->trackFitters().size();
  std::auto_ptr<TTTrackCollection> allFitTTTracksForOutput[nFitAlgs]; 
  map<string, unsigned int> locationInsideArray;
  unsigned int ialg = 0;
  for (const string& fitterName : settings_->trackFitters()) {
    std::auto_ptr<TTTrackCollection> fitTTTracksForOutput(new TTTrackCollection);
    allFitTTTracksForOutput[ialg] =  fitTTTracksForOutput;
    locationInsideArray[fitterName] = ialg++;
  }
  // Storage for EDM stub collections to be produced for comparison of hardware & software.
  std::auto_ptr<HwStubCollection>     outputSimStubs(new HwStubCollection);
  std::auto_ptr<HwStubCollection>  allOutputSimStubs(new HwStubCollection);
  std::auto_ptr<HwTrackCollection>         effTracks(new HwTrackCollection);
  std::auto_ptr<HwTrackCollection>     algoEffTracks(new HwTrackCollection);

  //=== Loop over matrix of Hough-Transform arrays, filling them with stubs.

  unsigned ntracks(0);
  // Fill Hough-Transform arrays with stubs.
  for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
    for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {

      Sector& sector = mSectors(iPhiSec, iEtaReg);
      HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);

      // Initialize constants for this sector.
      sector.init(settings_, iPhiSec, iEtaReg); 
      htPair.init(settings_, sector.etaMin(), sector.etaMax(), sector.phiCentre());

      for (const Stub* stub: vStubs) {
	// Digitize stub as would be at input to GP. This doesn't need the octant number, since we assumed an integer number of
	// phi digitisation  bins inside an octant. N.B. This changes the coordinates & bend stored in the stub.
	// The cast allows us to ignore the "const".
	if (settings_->enableDigitize()) (const_cast<Stub*>(stub))->digitizeForGPinput(iPhiSec);

	// Check if stub is inside this sector
        bool inside = sector.inside( stub );

        if (inside) {
	  // Check which eta subsectors within the sector the stub is compatible with (if subsectors being used).
	  const vector<bool> inEtaSubSecs =  sector.insideEtaSubSecs( stub );

	  // Digitize stub if as would be at input to HT, which slightly degrades its coord. & bend resolution, affecting the HT performance.
	  if (settings_->enableDigitize()) (const_cast<Stub*>(stub))->digitizeForHTinput(iPhiSec);

	  // Store stub in Hough transform array for this sector, indicating its compatibility with eta subsectors with sector.
	  htPair.store( stub, inEtaSubSecs );
	}
      }

      // Finish. Look for tracks in r-phi HT array etc.
      htPair.end();

      // Convert these tracks to EDM format for output (not used by Histos class).
      const vector<L1track3D>& vecTrk3D = htPair.trackCands3D();
      ntracks += vecTrk3D.size();
      for (const L1track3D& trk : vecTrk3D) {
        TTTrack< Ref_PixelDigi_ > htTTTrack = converter.makeTTTrack(trk, iPhiSec, iEtaReg);
        htTTTracksForOutput->push_back( htTTTrack );
      }
    }
  }

  // Initialize the duplicate track removal algorithm that can optionally be run after the track fit.
  KillDupFitTrks killDupFitTrks;
  killDupFitTrks.init(settings_, settings_->dupTrkAlgFit());
  
  //=== Do a helix fit to all the track candidates.

  vector<std::pair<std::string, L1fittedTrack>> fittedTracks;
  for (unsigned int iPhiSec = 0; iPhiSec < settings_->numPhiSectors(); iPhiSec++) {
    for (unsigned int iEtaReg = 0; iEtaReg < settings_->numEtaRegions(); iEtaReg++) {

      HTpair& htPair = mHtPairs(iPhiSec, iEtaReg);

      // In principal, should digitize stubs on track here using digitization relative to this phi sector.
      // However, previously digitized stubs will still be valid if digitization bins in phi align with
      // phi sector boundaries, so couldn't be bothered. If they do not align, the average effect of 
      // digitization on the track fit will be correct, but the effect on individual tracks will not.

      // Get track candidate sfound by Hough transform in this sector.
      const vector<L1track3D>& vecTrk3D = htPair.trackCands3D();
      // Loop over all the fitting algorithms we are trying.
      for (const string& fitterName : settings_->trackFitters()) {
        // Fit all tracks in this sector
	vector<L1fittedTrack> fittedTracksInSec;
        for (const L1track3D& trk : vecTrk3D) {
	  L1fittedTrack fitTrack = fitterWorkerMap_[fitterName]->fit(trk, iPhiSec, iEtaReg);
	  // Store fitted tracks, such that there is one fittedTracks corresponding to each HT tracks.
	  // N.B. Tracks rejected by the fit are also stored, but marked.
	  fittedTracksInSec.push_back(fitTrack);
	}

	// Run duplicate track removal on the fitted tracks if requested.
	// N.B. If a duplicate removal algorithm is run, it will also remove tracks rejected by the fitter.
	const vector<L1fittedTrack> filtFittedTracksInSec = killDupFitTrks.filter( fittedTracksInSec );

	// Store fitted tracks from entire tracker.
	for (const L1fittedTrack& fitTrk : filtFittedTracksInSec) {
	  fittedTracks.push_back(std::make_pair(fitterName, fitTrk));
	  // Convert these fitted tracks to EDM format for output (not used by Histos class).
	  // Only do this for valid fitted tracks, meaning that these EDM tracks do not correspond 1 to 1 with fittedTracks.
	  if (fitTrk.accepted()) {
	    TTTrack< Ref_PixelDigi_ > fitTTTrack = converter.makeTTTrack(fitTrk, iPhiSec, iEtaReg);
	    allFitTTTracksForOutput[locationInsideArray[fitterName]]->push_back(fitTTTrack);
	  }
	}
      }
    }
  }

  // Histogram the undigitized stubs, since with some firwmare versions, quantities like digitized stub dphi are
  // not available, so would give errors when people try histogramming them.
  for (const Stub* stub: vStubs) {
    if (settings_->enableDigitize()) (const_cast<Stub*>(stub))->reset_digitize();
  }

  //=== Fill histograms that check if choice of (eta,phi) sectors is good.
  hists_->fillEtaPhiSectors(inputData, mSectors);

  //=== Fill histograms that look at filling of r-phi HT arrays.
  hists_->fillRphiHT(mHtPairs);

  //=== Fill histograms that look at r-z filters (or other filters run after r-phi HT).
  hists_->fillRZfilters(mHtPairs);

  //=== Fill histograms studying track candidates found by r-phi Hough Transform.
  hists_->fillTrackCands(inputData, mSectors, mHtPairs);

  //=== Fill histograms studying track fitting performance
  hists_->fillTrackFitting(inputData, fittedTracks,  settings_->chi2OverNdfCut() );

  //=== Output digitized stubs in format expected by hardware for use by the comparison software,
  //=== which compares hardware with software.

  if (settings_->enableDigitize() && settings_->writeOutEdmFile()) {

    DemoOutput demoOutput(settings_);

    // Fill allOutputSimStubs and outputSimStubs with stubs stored in HardwareStub class.
    // The former contains all stubs; the latter only stubs assigned to L1 tracks.
    demoOutput.getStubCollection(mHtPairs,
   	   	                 allOutputSimStubs, outputSimStubs);

    // Fill effTracks and algoEffTracks with stubs on tracking particles stored in HardwareTrack class.
    // The former contains all TP, whilst the latter contains only those uses for algorithmic efficiency measurment.
    demoOutput.getTPstubCollection(mHtPairs, mSectors, vTPs,
				   effTracks, algoEffTracks);
  }

  //=== Store output EDM track and hardware stub collections.
  iEvent.put(htTTTracksForOutput,  "TML1TracksHT");
  for (const string& fitterName : settings_->trackFitters()) {
    string edmName = string("TML1Tracks") + fitterName;
    iEvent.put(allFitTTTracksForOutput[locationInsideArray[fitterName]], edmName);
  }
  iEvent.put(outputSimStubs,       "OutputSimStub");
  iEvent.put(allOutputSimStubs, "AllOutputSimStub");
  iEvent.put(effTracks,          "EfficiencyTrack");
  iEvent.put(algoEffTracks,  "AlgoEfficiencyTrack");
}


void TMTrackProducer::endJob() 
{
  hists_->endJobAnalysis();

  for (const string& fitterName : settings_->trackFitters()) {
    //cout << "# of duplicated stubs = " << fitterWorkerMap_[fitterName]->nDupStubs() << endl;
    delete fitterWorkerMap_[ string(fitterName) ];
  }

  cout<<endl<<"Number of (eta,phi) sectors used = (" << settings_->numEtaRegions() << "," << settings_->numPhiSectors()<<")"<<endl; 
}

DEFINE_FWK_MODULE(TMTrackProducer);
