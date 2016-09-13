#include "DataFormats/Demonstrator/interface/HardwareTrack.h"
#include <iostream>
#include <sstream>
#include <string>

l1t::HardwareTrack::HardwareTrack(unsigned int trackId, unsigned int phiSec,
      unsigned int etaReg, float pt, float eta, float phi0, int charge, int pdgId){
	trackId_ = trackId;
	phiSec_ = phiSec;
	etaReg_ = etaReg;
	pt_ = pt;
	eta_ = eta;
	phi0_ = phi0;
	charge_ = charge;
	pdgId_ = pdgId;
}

l1t::HardwareTrack::~HardwareTrack() 
{

}

 
