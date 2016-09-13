#ifndef DataFormats_L1Trigger_HardwareTrack_h
#define DataFormats_L1Trigger_HardwareTrack_h

#include <vector>
#include "DataFormats/Demonstrator/interface/HardwareStub.h"

namespace l1t {
  
  class HardwareStub;
  class HardwareTrack;
  /**
   * @brief     Class storing hardware tracks information at the output of the HT
   * @details   This class stores the information about track candidates found by the Hough Transform Fw and matched with the Monte Carlo simulation.
   * 
   * @author    Davide Cieri (davide.cieri@stfc.ac.uk)
   */

  class HardwareTrack{

  public:    
   /**
    * @brief      Normal HardwareTrack constructor
    *
    * @param[in]  trackId  The track identifier in the event
    * @param[in]  phiSec   The phi Sector in which the track has been found.
    * @param[in]  etaReg   The eta region in which the track has been found.
    * @param[in]  pt       The transverse momentum of the track
    * @param[in]  eta      The pseudo-rapidity coordinate of the track
    * @param[in]  phi0     The phi coordinate of the track at the production
    *                      point
    * @param[in]  charge   The charge
    * @param[in]  pdgId    The pdg id of the track
    */
    HardwareTrack(
      unsigned int trackId = 0,
      unsigned int phiSec=0,
      unsigned int etaReg=0,
      float pt=0,
      float eta = 0,
      float phi0=0,
      int charge = 0,
      int pdgId = 0
    );
    
    /// Normal Destructor
    ~HardwareTrack();

    /// The track identifier in the event
    unsigned int trackId()          const { return trackId_ ; }
    /// pdg Id of tracking particle
    int pdgId()                     const{ return pdgId_ ;       }
    /// Charge of the track
    int charge()                    const{ return charge_;       }
    /// Transverse momentum of the track
    float pt()                      const{ return pt_;           }
    /// Pseudo-rapidity of track
    float eta()                     const{ return eta_;          }
    /// Azimuthal angle phi of the track at the production point
    float phi0()                    const{ return phi0_; }
    /// Phi sector in which the track has been found.
    unsigned int phiSec()           const{ return phiSec_;       } 
    /// Eta Region in which the track has been found.
    unsigned int etaRegion()        const{ return etaReg_;      }
    
    /// Functions returning stubs associated with this track.
    const std::vector<l1t::HardwareStub>& assocStubs() const{ return assocStubs_; } 

    /// Number of associated stubs with the track
    unsigned int numAssocStubs()                    const{ return assocStubs_.size(); }
    
    /// Associate a stub to the track
    void Fill(l1t::HardwareStub stub) {assocStubs_.push_back(stub); }

  private:

    // Private members
    std::vector<l1t::HardwareStub> assocStubs_;
    unsigned int trackId_;
    unsigned int phiSec_;
    unsigned int etaReg_;
    float pt_;
    float eta_;
    float phi0_;
    int charge_;
    int pdgId_;

  };


  typedef std::vector<l1t::HardwareStub> HwTrackCollection;
}

#endif

