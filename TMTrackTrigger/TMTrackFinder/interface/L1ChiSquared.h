///=== This is the base class for the linearised chi-squared track fit algorithms.

///=== Written by: Sioni Summers and Alexander D. Morton

#ifndef __L1_CHI_SQUARED__
#define __L1_CHI_SQUARED__
 
#include "TMTrackTrigger/TMTrackFinder/interface/Matrix.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TrackFitGeneric.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"
#include <vector>
#include <map>
#include <utility>
 
 
class L1ChiSquared : public TrackFitGeneric{
public:
    L1ChiSquared(const Settings* settings, const uint nPar);
 
    virtual ~L1ChiSquared(){}
 
    L1fittedTrack fit(const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg);
 
    virtual std::string getParams()=0;
 
protected:
    /* Methods */
    virtual std::vector<double> seed(const L1track3D& l1track3D)=0;
    virtual std::vector<double> residuals(std::vector<double> x)=0;
    virtual Matrix<double> D(std::vector<double> x)=0; // derivatives
    virtual Matrix<double> Vinv()=0; // Covariances
    virtual std::map<std::string, double> convertParams(std::vector<double> x)=0;
 
    /* Variables */
    std::vector<const Stub*> stubs_;
    std::map<std::string, double> trackParams_;
    uint nPar_;
    float largestresid_;
    int ilargestresid_;
    double chiSq_;
 
private:

    void calculateChiSq( std::vector<double> resids );
    void calculateDeltaChiSq( std::vector<double> deltaX, std::vector<double> covX );

    int numFittingIterations_;
    int killTrackFitWorstHit_;
    double generalResidualCut_;
    double killingResidualCut_;

    unsigned int minStubLayers_;
    float minPtToReduceLayers_;
 
};
 
#endif
 

