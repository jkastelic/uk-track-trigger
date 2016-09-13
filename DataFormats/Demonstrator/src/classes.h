#include "DataFormats/Demonstrator/interface/HardwareStub.h"
#include "DataFormats/Demonstrator/interface/HardwareTrack.h"

#include "Rtypes.h" 
#include "Math/Cartesian3D.h" 
#include "Math/Polar3D.h" 
#include "Math/CylindricalEta3D.h" 
#include "Math/PxPyPzE4D.h" 
#include <boost/cstdint.hpp> 
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefProd.h"

namespace {
  struct dictionary {
    std::vector<l1t::HardwareStub> StubColl;
    std::vector<l1t::HwStubCollection> StubOverLinksColl;
    std::vector<l1t::HardwareTrack> HwTrackCollection;

    edm::Wrapper<std::vector<l1t::HardwareStub> > w_StubColl;
    edm::Wrapper<std::vector<l1t::HwStubCollection> > w_StubLinkColl;
    edm::Wrapper<std::vector<std::vector<l1t::HardwareStub> > > w_StubLinkColl0;

    edm::Wrapper<std::vector<l1t::HardwareTrack> > w_TrackColl;

   };
}
