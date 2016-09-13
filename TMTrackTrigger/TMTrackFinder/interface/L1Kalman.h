///=== This is the base class for the Kalman Filter track fit algorithm.

///=== Written by: Sioni Summers

#ifndef __L1_KALMAN__
#define __L1_KALMAN__
 
#include "TMTrackTrigger/TMTrackFinder/interface/Matrix.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TrackFitGeneric.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"
#include <map>
#include <vector>
 
class L1Kalman : public TrackFitGeneric{
 
    public:
        L1Kalman(const Settings* settings, const uint nPar);
 
        virtual ~L1Kalman(){}
 
        L1fittedTrack fit(const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg);
 
        virtual std::string getParams()=0;
 
    protected:
        /* Methods */
        virtual std::vector<double> seedx(const L1track3D& l1track3D)=0;
        virtual Matrix<double> seedP(const L1track3D& l1track3D)=0;
        virtual std::vector<double> d(const Stub* stub)=0;
        virtual Matrix<double> H(const Stub* stub)=0;
        virtual Matrix<double> F(const Stub* stub)=0;
        virtual Matrix<double> PxxModel()=0;
        virtual Matrix<double> PddMeas(const Stub* stub)=0;
        virtual std::map<std::string, double> convertParams(std::vector<double>)=0;
 
    private:
        unsigned nPar_;
 
};
#endif
 

