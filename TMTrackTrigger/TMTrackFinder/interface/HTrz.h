#ifndef __HTrz_H__
#define __HTrz_H__

#include "TMTrackTrigger/TMTrackFinder/interface/HTbase.h"

#include <vector>
#include <utility>

class Settings;
class Stub;
class TP;
class L1fittedTrack;


//=== The r-z Hough Transform array for a single (eta,phi) sector.
//===
//=== Its axes are (z0, zTrk), where zTrk is the z at which the track crosses a 
//=== user-configurable radius from the beam-line.


class HTrz : public HTbase {

public:
  
  HTrz() : HTbase() {}
  ~HTrz(){}

  // Initialization with cfg params, eta range covered by sector, and estimated q/Pt from previously run r-phi HT.
  void init(const Settings* settings, float etaMinSector, float etaMaxSector, float qOverPt);

  // Add stub to HT array.
  void store(const Stub* stub);

  // Termination. Causes HT array to search for tracks etc.
  // ... function end() is in base class ...

  //=== Info about track candidates found.

  // ... is available via base class ...

  //=== Utilities

  // Get the values of the track helix params corresponding to middle of a specified HT cell (i,j).
  // The helix parameters returned will be those corresponding to the two axes of the HT array.
  // So they might be (z0, zTrk), (z0, tan_lambda) or (z0, eta) etc. depending on the configuration.
  std::pair<float, float> helix2Dhough       (unsigned int i, unsigned int j) const;

  // Get the values of the track helix params corresponding to middle of a specified HT cell (i,j).
  // The helix parameters returned will be always be (z0, tan_lambda), irrespective of how the axes
  // of the HT array are defined.
  std::pair<float, float> helix2Dconventional(unsigned int i, unsigned int j) const;

  // Which cell in HT array should this TP be in, based on its true trajectory?
  // Returns (-1,-1) if TP not expected to be in any cell in this array.
  std::pair<int, int> trueCell( const TP* tp ) const;

  // Which cell in HT array should this fitted track be in, based on its fitted trajectory?
  // Returns (-1,-1) if fitted track not expected to be in any cell in this array.
  std::pair<int, int> getCell( const L1fittedTrack* fitTrk ) const;

  //--- Functions to check that stub filling is compatible with limitations of firmware.

  // Calculate maximum |gradient| that any stub's line across any of the r-phi HT arrays could have, to check it is < 1.
  static float maxLineGrad() {return maxLineGradient_;}
  // Summed over all r-phi HT arrays, returns fraction of stubs added to an HT column which do not lie NE, E or SE of stub added to previous HT column.
  static float fracErrorsTypeA() {return numErrorsTypeA_/float(numErrorsNormalisation_);}
  // Summed over all r-phi HT arrays, returns fraction of stubs added to more than 2 cells in one HT column. (Only a problem for Thomas' firmware).
  static float fracErrorsTypeB() {return numErrorsTypeB_/float(numErrorsNormalisation_);}

private:

  // For a given z0 bin, find the range of zTrk bins that a given stub is consistent with.
  std::pair<unsigned int, unsigned int> iZtrkRange( const Stub* stub, unsigned int iZ0Bin, bool debug = false) const;

  // Check that limitations of firmware would not prevent stub being stored correctly in this HT column.
  void countFirmwareErrors(unsigned int iZ0Bin, unsigned int iZtrkBinMin, unsigned int iZtrkBinMax);

  // Calculate maximum |gradient| that any stub's line across this HT array could have, so can check it doesn't exceed 1.
  float calcMaxLineGradArray() const;

  // Note if this is an r-phi or r-z Hough transform?
  bool isRphiHT() const {return false;}

  // Define the order in which the hardware processes rows of the HT array when it outputs track candidates.
  // Currently outputs an empty vector, meaning don't care.
  std::vector<unsigned int> rowOrder(unsigned int numRows) const {return std::vector<unsigned int>();}

private:

  //float trackerOuterRadius_; // Tracker outer radius.  

  // Specifications of HT array.

  unsigned int nBinsZ0Axis_;       // Number of bins in HT array in z0.
  float        maxAbsZ0Axis_;      // # Max. value of |z0| in HT array. Corresponds to beam-spot length.
  float        binSizeZ0Axis_;     // HT array bin size in z0.

  float        chosenRofZ_;      // Use z of track at radius="chosenRofZ" to define one of the r-z HT axes.
  unsigned int nBinsZtrkAxis_;   // Number of bins in HT array in zTrk
  float        minZtrkAxis_;     // Lower end of range in zTrk in HT array.
  float        maxZtrkAxis_;     // Upper end of range in zTrk in HT array.
  float        binSizeZtrkAxis_; // HT array bin size in zTrk.

  // Options when filling HT array.

  unsigned int killSomeHTCellsRz_; // Take all cells in HT array crossed by line corresponding to each stub (= 0) or take only some to reduce rate at cost of efficiency ( > 0)
  bool handleStripsRzHT_; // Should algorithm allow for uncertainty in stub (r,z) coordinate caused by length of 2S module strips when fill stubs in r-z HT?

  //--- Checks that stub filling is compatible with limitations of firmware.

  // Maximum |gradient| of line corresponding to any stub. Should be less than the value of 1.0 assumed by the firmware.
  static float maxLineGradient_;
  // Error count when stub added to cell which does not lie NE, E or SE of stub added to previous HT column.
  static unsigned int numErrorsTypeA_;
  // Error count when stub added to more than 2 cells in one HT column (problem only for Thomas' firmware).
  static unsigned int numErrorsTypeB_;
  // Error count normalisation
  static unsigned int numErrorsNormalisation_;

  // ... The Hough transform array data is in the base class ...

  // ... The list of found track candidates is in the base class ...
};
#endif
