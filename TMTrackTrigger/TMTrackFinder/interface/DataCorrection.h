#ifndef __DATACORRECTION_H__
#define __DATACORRECTION_H__

#include "FWCore/Utilities/interface/Exception.h"

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <cmath>



class DataCorrection {

  /*
   * Implementation of reduced bits for the bend information, 3 bits for PS, 4 bits for 2S.
   * This was originally proposed by Seb Viret in
   * https://espace.cern.ch/Tracker-Upgrade/Electronics/CIC/Shared%20Documents/Simulation%20studies/Bend_Codes_in_CIC_220415.pdf .
   * He did this encoding for the "2 GeV tight" stub window cut, designed to give 95% efficiency for
   * all stubs of Pt = 2 GeV. 
   *
   * However, this window size is no longer the recommended one by Seb Viret. The new recommended window size are narrower 
   * and are given in section "TTStubAlgorithm_tab2013_PixelDigi" of http://sviret.web.cern.ch/sviret/docs/CMS/SLHC/STUBS_base.py ,
   * The latter is an updated version of L1Trigger/TrackTrigger/python/TTStubAlgorithmRegister_cfi.py.
   * Seb Viret has not provided an update to his stub bend encoding that corresponds to these narrower windows, and recommends
   * that we do it ourselves, to profit from this narrower windows.
   * 
   * Ian Tomalin has therefore set the window sizes below according to those produced by Stub::frontendPass(),
   * used with a bend filter cut of 1.30. He noted the window sizes obtained using histogram
   * "hisBendFEVsLayerOrRing" produced by the "Histos" class. He then chose the encoding to fit within these
   * windows, with no serious attempt at tuning them. After implementing the encoding, he checked that
   * histogram "StubKillDataCorr" confirmed that few stubs were outside the windows assumed by the encoding.
   * 
   * IMPORTANT: If the code below detects a stub with bend outside the assumed window, meaning that looser windows
   * were used when generating the MC, then it sets boolean reject = true to tell you this. You should reject this stub.
   * You will effectively then be using tighter windows than originally used when producing the MC.
   */

public:

  //--- Given the original bend and barrel layer number,
  //--- return the degraded stub bend, a boolean indicatng if stub bend was outside the assumed window
  //--- size programmed below, and an integer indicating how many values of the original bend
  //--- were grouped together into this single value of the degraded bend.
  static void ConvertBarrelBend(float bend, unsigned int layer,
		float& degradedBend, bool& reject, unsigned int& num)
	{
		using namespace std;
		
    // Get degraded bend value
    DataCorrection::ConvertBarrelBendWork(bend, layer, degradedBend, reject);

    // Determine number of bend values that would lead to same degraded bend value.
    // This helps understand the loss in bend resolution caused by the bit encoding.

    // Stores number of bend values merged into a single one by bit encoding, for each layer & bend value.
    static std::vector< std::map<float, unsigned int> > barrelNums(30); // Dimension to larger than number of barrel laye
    std::map<float, unsigned int>& storedNum = barrelNums[layer];
		
    if (storedNum.find(degradedBend) != storedNum.end()) {

      // Calculation already done, so just retrieve it.
      num = storedNum[degradedBend];

    } else {

      // Calculation not yet done, so do it.

      num = 0;

      const int maxI = 30; // A number larger than the number of unique bend values.
      int maxAcceptedI = -1;
      float degradedBendI;
      bool rejectI;
      std::set<float> uniqueDegradedBends;
      for (int i = -maxI; i <= maxI; i++) {
	float bendI = 0.5*float(i);
	DataCorrection::ConvertBarrelBendWork(bendI, layer, degradedBendI, rejectI);
	if ( ! rejectI) {
	  if (abs(maxAcceptedI) < abs(i)) maxAcceptedI = abs(i);
	  if (degradedBend == degradedBendI) num++;
	  uniqueDegradedBends.insert(degradedBendI);
	}
      }

      //cout<<"Barrel layer="<<layer<<" Bend="<<bend<<" --> "<<degradedBend<<" reject=" <<reject<<" num="<<num<<endl;

      // Store result
      storedNum[degradedBend] = num;
    
      //--- Sanity checks
      if (maxAcceptedI < 0 || maxAcceptedI == (unsigned int) maxI) throw cms::Exception("DataCorrection:: barrel stub window size wrong. ")<<layer<<" "<<maxAcceptedI<<endl;
      // Number of degraded bend values should correspond to 3 bits (PS modules) or 4 bits (2S modules),
      // minus one, where the latter is because the encoding must be symmetric about 0.
      // Or perhaps less if no bit encoding was required.
      unsigned int numDegradedBendsExp = (layer <= 3)  ?  pow(2,3) - 1  :  pow(2,4) - 1;
      numDegradedBendsExp = min(numDegradedBendsExp, (unsigned int)(2*maxAcceptedI + 1)); 
      if (uniqueDegradedBends.size() != numDegradedBendsExp) throw cms::Exception("DataCorrection:: barrel stub encoding corresponds to wrong number of bits. ")<<layer<<" "<<numDegradedBendsExp<<" "<<uniqueDegradedBends.size()<<endl;
    }
  }

  //--- Given the original bend and endcap ring number,
  //--- return the degraded stub bend, a boolean indicating if stub bend was outside the assumed window
  //--- size programmed below, and an integer indicating how many values of the original bend
  //--- were grouped together into this single value of the degraded bend.

  static void ConvertEndcapBend(float bend, unsigned int ring,
   	float& degradedBend, bool& reject, unsigned int& num)
	{
    using namespace std;

    // Get degraded bend value
    DataCorrection::ConvertEndcapBendWork(bend, ring, degradedBend, reject);

    // Determine number of bend values that would lead to same degraded bend value.
    // This helps understand the loss in bend resolution caused by the bit encoding.

    // Stores number of bend values merged into a single one by bit encoding, for each ring & bend value
    static std::vector< std::map<float, unsigned int> > endcapNums(30); // Dimension to larger than number of endcap ring
    std::map<float, unsigned int>& storedNum = endcapNums[ring];

    if (storedNum.find(degradedBend) != storedNum.end()) {

      // Calculation already done, so just retrieve it.
      num = storedNum[degradedBend];

    } else {

      // Calculation not yet done, so do it.

      num = 0;

      const int maxI = 30; // A number larger than the number of unique bend values.
      int maxAcceptedI = -1;
      float degradedBendI;
      bool rejectI;
      set<float> uniqueDegradedBends;
      for (int i = -maxI; i <= maxI; i++) {
	float bendI = 0.5*float(i);
	DataCorrection::ConvertEndcapBendWork(bendI, ring, degradedBendI, rejectI);
	if ( ! rejectI) {
	  if (abs(maxAcceptedI) < abs(i)) maxAcceptedI = abs(i);
	  maxAcceptedI = i;
	  if (degradedBend == degradedBendI) num++;
	  uniqueDegradedBends.insert(degradedBendI);
	}
      }

      //cout<<"Endcap ring="<<ring<<" Bend="<<bend<<" --> "<<degradedBend<<" reject=" <<reject<<" num="<<num<<endl;

      // Store result
      storedNum[degradedBend] = num;
    
      //--- Sanity checks
      if (maxAcceptedI < 0 || maxAcceptedI == (unsigned int) maxI) throw cms::Exception("DataCorrection:: endcap stub window size wrong. ")<<ring<<" "<<maxAcceptedI<<endl;
      // Number of degraded bend values should correspond to 3 bits (PS modules) or 4 bits (2S modules),
      // minus one, where the latter is because the encoding must be symmetric about 0.
      // Or perhaps less if no bit encoding was required.
      unsigned int numDegradedBendsExp = (ring <= 9)  ?  pow(2,3) - 1  :  pow(2,4) - 1;
      numDegradedBendsExp = min(numDegradedBendsExp, (unsigned int)(2*maxAcceptedI + 1)); 
      if (uniqueDegradedBends.size() != numDegradedBendsExp) throw cms::Exception("DataCorrection:: endcap stub encoding corresponds to wrong number of bits. ")<<ring<<" "<<numDegradedBendsExp<<" "<<uniqueDegradedBends.size()<<endl;
    }
  }

private:

  // No constructor needed, since all function members are static.
  DataCorrection() = delete;

  //--- Given the original bend and barrel layer number,
  //--- return the degraded stub bend & a boolean indicating if stub bend was outside the assumed window
  //--- size programmed below.

  static void ConvertBarrelBendWork(float bend, unsigned int layer,
		float& degradedBend, bool& reject)
	{
    using namespace std;
    
    float b = fabs(bend);
    int   s = (bend > 0) ? 1 : -1;

    degradedBend = b; // default to non-degraded.
    reject = true;

    // Layers 1-3 have PS modules & layers 4-6 have 2S modules
    switch( layer ){
    case 1:
    case 2:
	if (b <= 1.5) reject = false; 
        break;
    case 3:
	if (b <= 2.5) reject = false; 
        if      (1.0 <= b && b <= 1.5) degradedBend = 1.25;
        else if (2.0 <= b && b <= 2.5) degradedBend = 2.25;
        break;
    case 4:
	if (b <= 3.5) reject = false; 
        break;
    case 5:
	if (b <= 4.5) reject = false; 
        if      (3.0 <= b && b <= 3.5) degradedBend = 3.25;
        else if (4.0 <= b && b <= 4.5) degradedBend = 4.25;
        break;
    case 6:
	if (b <= 5.0) reject = false; 
        if      (2.5 <= b && b <= 3.0) degradedBend = 2.75;
        else if (3.5 <= b && b <= 4.0) degradedBend = 3.75;
        else if (4.5 <= b && b <= 5.0) degradedBend = 4.75;
        break;
    default:
        throw cms::Exception("DataCorrection:: unknown barrel layer.")<<layer<<endl;
        break;
    }

    if (reject) degradedBend = 99999.;

    degradedBend *= s;
  }

  //--- Given the original bend and endcap ring number,
  //--- return the degraded stub bend & a boolean indicating if stub bend was outside the assumed window
  //--- size programmed below.

  static void ConvertEndcapBendWork(float bend, unsigned int ring,
		float& degradedBend, bool& reject)
	{
    using namespace std;
    
    float b = fabs(bend);
    int   s = (bend > 0) ? 1 : -1;

    degradedBend = b; // default to non-degraded.
    reject = true;

    // Rings 1-9 have PS modules & layers 10-15 have 2S modules
    switch( ring ){
    case 1:
    case 2:
	if (b <= 1.0) reject = false; 
        break;
    case 3:
	if (b <= 1.5) reject = false; 
        break;
    case 4:
    case 5:
	if (b <= 2.0) reject = false; 
        if      (1.5 <= b && b <= 2.0) degradedBend = 1.75;
        break;
    case 6:
    case 7:
	if (b <= 2.5) reject = false; 
        if      (1.0 <= b && b <= 1.5) degradedBend = 1.25;
        else if (2.0 <= b && b <= 2.5) degradedBend = 2.25;
        break;
    case 8:
    case 9:
	if (b <= 3.0) reject = false; 
        if      (0.5 <= b && b <= 1.0) degradedBend = 0.75;
        else if (1.5 <= b && b <= 2.0) degradedBend = 1.75;
        else if (2.5 <= b && b <= 3.0) degradedBend = 2.75;
        break;
    case 10:
	if (b <= 4.0) reject = false; 
        if      (3.5 <= b && b <= 4.0) degradedBend = 3.75;
        break;
    case 11:
    case 12:
    case 13:
	if (b <= 3.5) reject = false; 
        break;
    case 14:
	if (b <= 4.0) reject = false; 
        if      (3.5 <= b && b <= 4.0) degradedBend = 3.75;
        break;
    case 15:
	if (b <= 4.5) reject = false; 
        if      (3.0 <= b && b <= 3.5) degradedBend = 3.25;
        else if (4.0 <= b && b <= 4.5) degradedBend = 4.25;
        break;
    default:
        throw cms::Exception("DataCorrection:: unknown endcap ring.")<<ring<<endl;
        break;
    }

    if (reject) degradedBend = 99999.;

    degradedBend *= s;
  }
};

#endif

