///=== This is the linearized chi2 track fit algorithm.
///=== It is based on the tracklet group's track fit software 
///=== https://github.com/EmanuelPerez/cmssw/blob/TTI_62X_TrackTriggerObjects/SimDataFormats/SLHC/interface/L1TTrack.hh

///=== Written by: Alexander D. Morton


#ifndef __TrackFitLinearAlgo_H__
#define __TrackFitLinearAlgo_H__

// Don't fit track candidates if they have more than this number of stubs.
#define __MAX_STUBS_PER_TRK__   30

#include "TMTrackTrigger/TMTrackFinder/interface/TrackFitGeneric.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"

#include <vector>
#include <utility>


class Settings;

class TrackFitLinearAlgo: public TrackFitGeneric {

public:

  // Set configuration parameters.
 TrackFitLinearAlgo( const Settings* settings, const uint nPar );

  ~TrackFitLinearAlgo() {}

  // Initialize contants at start of run.
  void initRun();

  // Fit a track candidate obtained from the Hough Transform.
  // Specify which phi sector and eta region it is in.
   L1fittedTrack fit( const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg );

   std::string getParams();

private:
  // Method to invert matrix by hand
  void invert( float M[5][10], unsigned int n );

  // Method to calculate derivatives
  void calculateDerivatives( bool withd0 );

  // Method to calculate residuals
  void residuals( float& largestresid, int& ilargestresid );

  // Linear track fitting method
  void linearTrackFit( bool withd0 );

private:
  // Configuration parameters
  const Settings* settings_;
  float invPtToInvR_;
  uint nPar_;

  const bool print;

  std::vector< const Stub* > stubs_;
  float rinv_;
  float phi0_;
  float z0_;
  float t_;
  float d0_;

  float rinvfit4par_;
  float phi0fit4par_;
  float z0fit4par_;
  float tfit4par_;

  float rinvfit_;
  float phi0fit_;
  float d0fit_;
  float z0fit_;
  float tfit_;

  float chisq_;
  float chisq4par_;

  int ichisq1_;
  int ichisq2_;

  float D_[5][2*__MAX_STUBS_PER_TRK__];
  
  float M_[5][10];
  
  float MinvDt_[5][2*__MAX_STUBS_PER_TRK__];

  // Configuration parameters
  int numFittingIterations_;
  bool killTrackFitWorstHit_;
  double generalResidualCut_;
  double killingResidualCut_;
  unsigned int minStubLayers_;
  float minPtToReduceLayers_;

  std::string configParameters_;

};
#endif

