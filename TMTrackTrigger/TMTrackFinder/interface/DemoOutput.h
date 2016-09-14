#ifndef __DEMOOUTPUT_H__
#define __DEMOOUTPUT_H__

#include <DataFormats/Demonstrator/interface/HardwareStub.h>
#include <DataFormats/Demonstrator/interface/HardwareTrack.h>
#include <TMTrackTrigger/TMTrackFinder/interface/HTpair.h>
#include <TMTrackTrigger/TMTrackFinder/interface/Sector.h>
#include <TMTrackTrigger/TMTrackFinder/interface/TP.h>

#include <memory>

class Settings;


//==================================================================================================
/**
* Make EDM stub collections (in "HardwareStub" class) needed for comparison of hardware & software.
* Davide Cieri is contact person for contents of this class.
*/
//==================================================================================================

class DemoOutput {

public:

	typedef std::vector<l1t::HardwareStub>           HwStubCollection;
	typedef std::vector<l1t::HardwareTrack>          HwTrackCollection;

	DemoOutput(const Settings* settings) : settings_(settings) {}

	// Create EDM collections of all stubs and of the subset of stubs that are on L1 tracks, stored in the HardwareStub class.
	void getStubCollection(
		const boost::numeric::ublas::matrix<HTpair>& mHtPairs,
		std::auto_ptr<HwStubCollection>& hwStubs,
		std::auto_ptr<HwStubCollection>& hwStubsOnTracks
	);

	// Create EDM collectiosn of all stubs produced by all tracking particles, or by the subset of tracking particles used
	// for the algorithmic tracking efficiency measurement. Store the results in the HardwareTrack class.
	void getTPstubCollection(
		const boost::numeric::ublas::matrix<HTpair>& mHtPairs,
		const boost::numeric::ublas::matrix<Sector>& mSectors,
		const std::vector<TP>& vTPs,
		std::auto_ptr<HwTrackCollection>& hwStubsOnTP,
		std::auto_ptr<HwTrackCollection>& hwStubsOnTPforAlgEff
	);

	~DemoOutput() {}

private:

	const Settings* settings_;
};
#endif

