#include "DataFormats/Demonstrator/interface/HardwareStub.h"
#include <iostream>
#include <sstream>
#include <string>

using namespace l1t;

// Normal empty constructor
HardwareStub::HardwareStub() 
{
  phiSec_ = 0;
  eta_reg_ = 0;
  dummy_ = 0;
  valid_ = 0;
  rho_ = 0;
  m_bin_ = 0;
  c_bin_ = 0;
  phiS_ = 0;
  rT_ = 0;
  z_ = 0;
  dphi_ = 0;
  mbin_min_ = 0;
  mbin_max_ = 0;
  layerId_ = 0;
  phiOctant_ = 0;
  phiO_ = 0;
  bend_ = 0;
  moduleType_ = 0;
}

// Constructor for the new HT fw input DataFormat (DaisyChain)
HardwareStub::HardwareStub( unsigned int phiSec,             
      unsigned int eta_reg,       
      bool dummy,
      bool valid,
      int phiS,                   
      int rT,                     
      int z,                      
      int mMin,                   
      int mMax,                   
      unsigned int layerId,
      unsigned int subeta,
      bool deadLayer,
      unsigned int deadLayerId){
  phiSec_ = phiSec;
  eta_reg_ = eta_reg;
  dummy_ = dummy;
  valid_ = valid;
  phiS_ = phiS;
  rT_ = rT;
  z_ = z;
  mbin_min_ = mMin;
  mbin_max_ = mMax;
  layerId_ = layerId;
  subeta_ = subeta;
  deadLayer_ = deadLayer;
  deadLayerId_ = deadLayerId_;
}

// Constructor for the new HT fw output DataFormat (DaisyChain)
HardwareStub::HardwareStub(
      unsigned int phiSec,             
      unsigned int eta_reg,       
      bool dummy,
      bool valid,
      unsigned int layerId,
      int phiS,                   
      int rT,                     
      int z,                      
      int cBin,                   
      int mBin              
  ){
  phiSec_ = phiSec;
  eta_reg_ = eta_reg;
  dummy_ = dummy;
  valid_ = valid;
  layerId_ = layerId;
  phiS_ = phiS;
  rT_ = rT;
  z_ = z;
  c_bin_ = cBin;
  m_bin_ = mBin;
}

// Constructor for Old HT input DataFormats
HardwareStub::HardwareStub(
    unsigned int phiSec,             
    unsigned int eta_reg,       
    bool dummy,
    bool valid,
    unsigned int rho,
    int phiS,                   
    int rT,                     
    int z,                      
    int dphi
  ){
  phiSec_ = phiSec;
  eta_reg_ = eta_reg;
  dummy_ = dummy;
  valid_ = valid;
  rho_ = rho;
  phiS_ = phiS;
  rT_ = rT;
  z_ = z;
  dphi_ = dphi;
}


// Constructor for full working GP input DataFormat
HardwareStub::HardwareStub( unsigned int phiOctant, int phiO, int rT, int z, int bend, unsigned int moduleType, bool valid, bool dummy, unsigned int stubId){
  phiOctant_ = phiOctant;
  dummy_ = dummy;
  valid_ = valid;
  phiO_ = phiO;
  rT_ = rT;
  z_ = z;
  bend_ = bend;
  moduleType_ = moduleType;
  stubId_ = stubId;
  phiSec_ = 0;
}

// Constructor for half working GP input DataFormat
HardwareStub::HardwareStub( unsigned int phiOctant, int phiO, int rT, std::vector<bool> vEtaRegions, std::vector<bool> vPhiSectors, unsigned int layerId, int mMin, int mMax, bool valid, bool dummy, unsigned int stubId ){
  
  phiOctant_ = phiOctant;
  dummy_ = dummy;
  valid_ = valid;
  phiO_ = phiO;
  rT_ = rT;
  vEtaRegions_ = vEtaRegions;
  vPhiSectors_ = vPhiSectors;
  mbin_min_ = mMin;
  mbin_max_ = mMax;
  stubId_ = stubId;
  layerId_ = layerId;
}

// Destructor 
HardwareStub::~HardwareStub() 
{

}

 
