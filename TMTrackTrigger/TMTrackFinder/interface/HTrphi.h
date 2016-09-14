#ifndef __HTrphi_H__
#define __HTrphi_H__

#include "TMTrackTrigger/TMTrackFinder/interface/HTbase.h"

#include <vector>
#include <utility>

class Settings;
class Stub;
class TP;
class L1fittedTrack;


//=== The r-phi Hough Transform array for a single (eta,phi) sector.
//===
//=== Its axes are (q/Pt, phiTrk), where phiTrk is the phi at which the track crosses a 
//=== user-configurable radius from the beam-line.

class HTrphi : public HTbase {

public:
  
  HTrphi() : HTbase() {}
  ~HTrphi(){}

  // Initialization with eta range covered by sector and phi coordinate of its centre.
  void init(const Settings* settings, float etaMinSector, float etaMaxSector, float phiCentreSector);

  // Add stub to HT array.
  // If eta subsectors are being used within each sector, specify which ones the stub is compatible with.
  void store( const Stub* stub, const std::vector<bool>& inEtaSubSecs);

  // Termination. Causes HT array to search for tracks etc.
  // ... function end() is in base class ...

  //=== Info about track candidates found.

  // ... is available via base class ...

  //=== Utilities

  // Get the values of the track helix params corresponding to middle of a specified HT cell (i,j).
  // The helix parameters returned will be those corresponding to the two axes of the HT array.
  // So they might be (q/pt, phi0) or (q/pt, phi65) etc. depending on the configuration.
  
  std::pair<float, float> helix2Dhough       (unsigned int i, unsigned int j) const;

  // Get the values of the track helix params corresponding to middle of a specified HT cell (i,j).
  // The helix parameters returned will be always be (q/pt, phi0), irrespective of how the axes
  // of the HT array are defined.
  std::pair<float, float> helix2Dconventional(unsigned int i, unsigned int j) const;

  // Which cell in HT array should this TP be in, based on its true trajectory?
  // Returns (-1,-1) if TP not expected to be in any cell in this array.
  std::pair<int, int> trueCell( const TP* tp ) const;

  // Which cell in HT array should this fitted track be in, based on its fitted trajectory?
  // Returns (-1,-1) if fitted track not expected to be in any cell in this array.
  std::pair<int, int> getCell( const L1fittedTrack* fitTrk ) const;

  // Check if specified cell has been merged with its 2x2 neighbours into a single cell,
  // as it is in low Pt region.
  bool mergedCell(unsigned int iQoverPtBin, unsigned int jPhiTrkBin) const;

  //--- Functions to check that stub filling is compatible with limitations of firmware.

  // Calculate maximum |gradient| that any stub's line across any of the r-phi HT arrays could have, to check it is < 1.
  static float maxLineGrad() {return maxLineGradient_;}
  // Summed over all r-phi HT arrays, returns fraction of stubs added to an HT column which do not lie NE, E or SE of stub added to previous HT column.
  static float fracErrorsTypeA() {return numErrorsTypeA_/float(numErrorsNormalisation_);}
  // Summed over all r-phi HT arrays, returns fraction of stubs added to more than 2 cells in one HT column. (Only a problem for Thomas' firmware).
  static float fracErrorsTypeB() {return numErrorsTypeB_/float(numErrorsNormalisation_);}

private:

  // For a given Q/Pt bin, find the range of phi bins that a given stub is consistent with.
  std::pair<unsigned int, unsigned int> iPhiRange( const Stub* stub, unsigned int iQoverPtBin, bool debug = false) const;

  // Check that limitations of firmware would not prevent stub being stored correctly in this HT column.
  void countFirmwareErrors(unsigned int iQoverPtBin, unsigned int iPhiTrkBinMin, unsigned int iPhiTrkBinMax);

  // Calculate maximum |gradient| that any stub's line across this HT array could have, so can check it doesn't exceed 1.
  float calcMaxLineGradArray() const;

  // Note if this is an r-phi or r-z Hough transform?
  bool isRphiHT() const {return true;}

  // Define the order in which the hardware processes rows of the HT array when it outputs track candidates.
  // Currently corresponds to highest Pt tracks first.
  // If two tracks have the same Pt, the -ve charge one is output before the +ve charge one.
  std::vector<unsigned int> rowOrder(unsigned int numRows) const;

private:

  float invPtToDphi_; // conversion constant.

  //--- Specifications of HT array.

  float maxAbsQoverPtAxis_;       // Max. |q/Pt| covered by  HT array.
  unsigned int nBinsQoverPtAxis_; // Number of bins in HT array in q/Pt.
  float binSizeQoverPtAxis_;      // HT array bin size in q/Pt.

  float chosenRofPhi_;             // Use phi of track at radius="chosenRofPhi" to define one of the r-phi HT axes.
  float phiCentreSector_;         // phiTrk angle of centre of this (eta,phi) sector.
  float maxAbsPhiTrkAxis_;        // Half-width of phiTrk axis in HT array.
  unsigned int nBinsPhiTrkAxis_;  // Number of bins in HT array in phiTrk axis.
  float binSizePhiTrkAxis_;       // HT array bin size in phiTrk

  bool  enableMerge2x2_; // Optionally merge 2x2 neighbouring cells into a single cell at low Pt, to reduce efficiency loss due to scattering.
  float minInvPtToMerge2x2_;

  //--- Options when filling HT array.

  unsigned int killSomeHTCellsRphi_; // Take all cells in HT array crossed by line corresponding to each stub (= 0) or take only some to reduce rate at cost of efficiency ( > 0)
  bool handleStripsRphiHT_; // Should algorithm allow for uncertainty in stub (r,z) coordinate caused by length of 2S module strips when fill stubs in r-phi HT?

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

