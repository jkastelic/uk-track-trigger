#ifndef DataFormats_L1Trigger_HardwareStub_h
#define DataFormats_L1Trigger_HardwareStub_h

#include <vector>

namespace l1t {
  class Cell;
  class HardwareStub;
  /**
   * @brief     Class for hardware stub at the different stages of the TMT hardware demonstrator
   * @details   This class is the representation in hardware of the stub coordinates. Different constructors are available for the different stage of the demonstrator. For an overview of the utilised DataFormats have a look at https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/demonstrator/specifications/demonstrator_formats_working_doc.docx
   * 
   * 
   * @author    Davide Cieri (davide.cieri@stfc.ac.uk)
   */
  class HardwareStub{

  public:    
    
    /**
     * @brief      Empty Constructor
     */
    HardwareStub();

    /**
     * @brief      Normal Constructor for the hardware stub with new HT input
     *             DataFormats (DaisyChain fw)
     *
     * @param[in]  phiSec       The phi Sector
     * @param[in]  eta_reg      The eta Region
     * @param[in]  dummy        The dummy bit
     * @param[in]  valid        The valid bit
     * @param[in]  phiS         The phiS (azimuthal angle phi of the stubs
     *                          relative to the phi sector)
     * @param[in]  rT           The rT value (radial position shifted by T)
     * @param[in]  z            z position
     * @param[in]  mMin         The lowest q/pT bin in which the stubs is
     *                          allowed to enter
     * @param[in]  mMax         The top q/pT bin in which the stubs is allowed
     *                          to enter
     * @param[in]  layerId      The tracker layer identifier
     * @param[in]  subeta       The Eta SubRegion identifier (default = 0 )
     * @param[in]  deadLayer    The dead layer bit
     * @param[in]  deadLayerId  The dead layer identifier
     */
    HardwareStub(
      unsigned int phiSec,             
      unsigned int eta_reg,       
      bool dummy,
      bool valid,
      int phiS,                   
      int rT,                     
      int z,                      
      int mMin,                   
      int mMax,                   
      unsigned int layerId,
      unsigned int subeta = 0,
      bool deadLayer = 0,
      unsigned int deadLayerId = 0
    );

    /**
     * @brief      Normal Constructor for the hardware stub with new HT output DataFormats (DaisyChain fw)
     *
     * @param[in]  phiSec   The phi Sector
     * @param[in]  eta_reg  The eta Region
     * @param[in]  dummy    The dummy bit
     * @param[in]  valid    The valid bit
     * @param[in]  layerId  The tracker layer identifier
     * @param[in]  phiS     The phiS (azimuthal angle phi of the stubs relative
     *                      to the phi sector)
     * @param[in]  rT       The rT value (radial position shifted by T)
     * @param[in]  z        z position
     * @param[in]  cBin     The c bin in the array (phi58 of the candidate track)
     * @param[in]  mBin     The m bin in the array (q/pT of the candidate track)
     */
    HardwareStub(
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
    );

    /**
     * @brief      Normal Constructor for the hardware stub with old HT input DataFormats (pipelined HT and systolic array firmware)
     *
     * @param[in]  phiSec   The phi Sector
     * @param[in]  eta_reg  The eta Region
     * @param[in]  dummy    The dummy bit
     * @param[in]  valid    The valid bit
     * @param[in]  rho      The rho value (bend filter precision)
     * @param[in]  phiS     The phiS (azimuthal angle phi of the stubs relative to the phi sector)
     * @param[in]  rT       The rT value (radial position shifted by T)
     * @param[in]  z        z position
     * @param[in]  dphi     Delta(phi) = bend angle of stub
     * @param[in]  layerId  The tracker layer identifier
     */
    HardwareStub(
      unsigned int phiSec,             
      unsigned int eta_reg,       
      bool dummy,
      bool valid,
      unsigned int rho,
      int phiS,                   
      int rT,                     
      int z,                      
      int dphi
    );

    /**
    * @brief Constructor for GP input DataFormats
    *
    * @param[in]  phiOctant   The phi octant
    * @param[in]  phiO        The phi coordinate of the stub relative to the centre of the octant
    * @param[in]  rT          The r coordinate of the stub shifted by T
    * @param[in]  z           z coordinate of the stub
    * @param[in]  bend        The bend
    * @param[in]  moduleType  The module type
    * @param[in]  valid       The valid bit
    * @param[in]  dummy       The dummy bit
    * @param[in]  stubId      The stub identifier
    */

    HardwareStub(
      unsigned int phiOctant,             
      int phiO,
      int rT,
      int z,
      int bend,
      unsigned int moduleType,
      bool valid,
      bool dummy,
      unsigned int stubId
    );

    /// Constructor for temp GP DataFormats
    HardwareStub( unsigned int phiOctant, int phiO, int rT, std::vector<bool> vEtaRegions, std::vector<bool> vPhiSectors, unsigned int layerId, int mMin, int mMax, bool valid, bool EoF, unsigned int stubId );    
  
    ~HardwareStub();

    /// Returns the phi Sector id
    unsigned int phiSec()           const{ return phiSec_;       }
    /// Returns the Eta Region id
    unsigned int etaRegion()        const{ return eta_reg_;      } 
    /// phi coordinate relative to the sector
    int phiS()                      const{ return phiS_;         } 
    /// r coordinate relative to the chosen radius
    int rT()                        const{ return rT_;           } 
    /// z coordinate
    int z()                         const{ return z_;            } 
    /// Delta(phi) = bend angle of stub
    int dphi()                      const{ return dphi_;         } 
    /// rho = Delta(phi) resolution
    unsigned int rho()              const{ return rho_;          }
     /// Reduced (3 bits) Delta(phi) coordinate
    int dphi_reduced()              const{ return dphi_reduced_; }
    /// q/pT bin of the candidate track 
    int m_bin()                     const{ return m_bin_;        }
    /// phiT bin of the candidate track relative to the sector 
    int c_bin()                     const{ return c_bin_;        }
    /// Bit declaring the dummy stub
    bool dummyStub()                const{ return dummy_;      }
     /// Bit declaring that the stub is valid
    bool valid()                    const{ return valid_;        }
    /// Tracker LayerId (reduced to 3 bits)
    unsigned int layerId()          const{ return layerId_;      }
    /// Lower q/pT bin in which the stub is allowed
    int MbinMin()                   const{ return mbin_min_;     }
    /// Upper q/pT bin in which the stub is allowed
    int MbinMax()                   const{ return mbin_max_;     }
    /// Debug Eta Region Number for GP fw
    unsigned int RealEtaRegion()    const{ return eta_reg_real_; }
    /// Debug PhiT Sector Number for GP fw
    unsigned int RealS()            const{ return S_real_;       }
    /// Bend
    int Bend()                      const{ return bend_;         }
    /// ModuleType ( 0 = PS_1 barrel, 1 = PS_2 barrel,  2 = PS endcap , 3 = 2S barrel 4 = endcap out , 5=  2s endcap in) 
    unsigned int ModuleType()       const{ return moduleType_;   }
    /// PhiOctant iD
    unsigned int phiOctant()        const{ return phiOctant_;    }
    /// Phi coordinate relative to the centre of the phi octant
    int phiO()                      const{ return phiO_;         } 
    /// Phi sector Id in the Octant (0-3)
    int SectorInOctant()            const{ return phiSec_ - phiOctant_*4; }
    /// Sub Eta sector (2 left, 1 right, 3 both, 0 null)
    unsigned int subEta()           const{ return subeta_;        }
    /// For GP debugging (Eta regions where the stub is shared )
    std::vector<bool> vEtaRegions()  const { return vEtaRegions_; }
    /// For GP debugging (phi sectors where the stub is shared )
    std::vector<bool> vPhiSectors()  const { return vPhiSectors_; }
    /// If true stub belongs to a dead layer in the tracker
    bool DeadLayer()                const { return deadLayer_;  }
    /// Id of the dead layer in the tracker (only if DeadLayer() == True)
    unsigned int DeadLayerId()      const { return deadLayerId_;}    
    // DEBUG PARAMETERS (FLOATING POINT)
    float fPhiS()                 const{ return fPhiS_; } 
    float fPhi()                  const{ return fPhi_;  }
    float fRt()                   const{ return fRt_;   }
    float fPhiO()                 const{ return fPhiO_; }
    float fZ()                    const{ return fZ_;    }
    unsigned int stubId()         const{ return stubId_;}
    bool RightPhiSec()      const{ return rightPhi_;}
    bool RightEtaReg()      const{ return rightEta_;}

    // Functions to set individually parameters inside the class
    void setPhiSec(unsigned int et)      { phiSec_ = et;        }
    void setEtaRegion(unsigned int et)   { eta_reg_ = et;       }
    void setphiS( int et)                { phiS_ = et;          }
    void setrT( int et )                 { rT_ = et;            }
    void setz( int et )                  { z_ = et;             }
    void setdphi( int et )               { dphi_ = et;          }
    void setrho(unsigned int et)         { rho_ = et;           }
    void setdphi_reduced(int et)         { dphi_reduced_ = et;  }
    void setm_bin(int et)                { m_bin_ = et;         }
    void setc_bin(int et)                { c_bin_ = et;         }
    void setDummyBit(unsigned int et)    { dummy_ = et;       }
    void setValidBit( unsigned int et )  { valid_ = et;         }
    void setLayerId(unsigned int et)     { layerId_ = et;       }
    void setMbinMin(int et)              { mbin_min_ = et;      }
    void setMbinMax(int et)              { mbin_max_ = et;      }
    void setRealEta(unsigned int et)     { eta_reg_real_ = et;  }
    void setRealS( unsigned int et )     { S_real_ = et;        }
    void setBend( int et )               { bend_ = et;          }
    void setPhiOctant( unsigned int et ) { phiOctant_ = et;     }
    void setModuleType( unsigned int et ){ moduleType_ = et;    }
    void setPhiO( int et )               { phiO_ = et;          }
    void setFphiS( float et )            { fPhiS_ = et;         }
    void setFphi( float et  )            { fPhi_ = et;          }
    void setFrT( float et   )            { fRt_ = et;           }
    void setFphiO(float et)              { fPhiO_ = et;         }
    void setFz(float et)                 { fZ_ = et;            }
    void setStubId(unsigned int et)      { stubId_ = et;        }
    void riseEtaRegBit( unsigned int et ){ vEtaRegions_[et] = 1; }
    void risePhiSecBit( unsigned int et ){ vPhiSectors_[et] = 1; }
    void PhiSecIsRight(bool et)          { rightPhi_ = et;       }
    void EtaRegIsRight(bool et)          { rightEta_ = et;       }
    void setVecEtaRegions( std::vector<bool> et) { vEtaRegions_ = et; }
    void setVecPhiSectors( std::vector<bool> et) { vPhiSectors_ = et; }

  private:
    // Private members
    unsigned int phiOctant_;
    unsigned int phiSec_;
    unsigned int eta_reg_;
    int phiS_;
    int phiO_;
    unsigned int subeta_;
    int rT_;
    int z_; 
    int dphi_;
    unsigned int rho_;
    int dphi_reduced_;
    int mbin_min_;
    int mbin_max_;
    unsigned int layerId_;
    int m_bin_;
    int c_bin_;
    int bend_;
    unsigned int moduleType_;
    unsigned int dummy_;
    unsigned int valid_;
    unsigned int eta_reg_real_;
    unsigned int S_real_;
    float fPhiS_;
    float fPhi_;
    float fRt_;
    float fPhiO_;
    float fZ_;
    unsigned int stubId_;
    std::vector<bool> vEtaRegions_;
    std::vector<bool> vPhiSectors_;
    bool rightEta_;
    bool rightPhi_;
    bool deadLayer_;
    unsigned int deadLayerId_;
  };

  typedef std::vector<HardwareStub> HwStubCollection;
  typedef std::vector<HwStubCollection> HwStubOverLinksCollection;
}

#endif

