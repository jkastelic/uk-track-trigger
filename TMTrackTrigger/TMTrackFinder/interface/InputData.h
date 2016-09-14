#ifndef __INPUTDATA_H__
#define __INPUTDATA_H__

#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include <vector>

class Settings;


//=== Unpacks stub & tracking particle (truth) data into user-friendlier format in Stub & TP classes.
//=== Also makes B-field available to Settings class.

class InputData
{
public:
	InputData(const edm::Event& iEvent, const edm::EventSetup& iSetup, Settings* settings);

	// Get tracking particles
	const std::vector<TP>&          getTPs()      const {return vTPs_;}
	// Get stubs that would be output by the front-end readout electronics 
	const std::vector<const Stub*>& getStubs()    const {return vStubs_;}

	//--- of minor importance ...

	// Get number of stubs prior to applying tighted front-end readout electronics cuts specified in section StubCuts of Analyze_Defaults_cfi.py. (Only used to measure the efficiency of these cuts).
	const std::vector<Stub>&        getAllStubs() const {return vAllStubs_;}

private:
	// tracking particles
	std::vector<TP> vTPs_;
	
	// stubs that would be output by the front-end readout electronics.
	std::vector<const Stub*> vStubs_;

	//--- of minor importance ...

	// all stubs, even those that would fail any tightened front-end readout 
	// electronic cuts specified in section StubCuts of Analyze_Defaults_cfi.py. 
	// (Only used to measure the efficiency of these cuts).
	std::vector<Stub> vAllStubs_;
};
#endif

