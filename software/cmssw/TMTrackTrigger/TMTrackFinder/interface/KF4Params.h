///=== This is the Kalman Filter for 4 helix parameters track fit algorithm.

///=== Written by: Sioni Summers

#ifndef __KF4PARAMS__
#define __KF4PARAMS__
 
#include "TMTrackTrigger/TMTrackFinder/interface/L1Kalman.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Matrix.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"
#include <map>
 
class KF4Params : public L1Kalman{
 
    public:
        KF4Params(const Settings* settings, const uint nPar);
        ~KF4Params(){}
        std::string getParams();
 
    protected:
        std::vector<double> seedx(const L1track3D& l1track3D);
        Matrix<double> seedP(const L1track3D& l1track3D);
        std::vector<double> d(const Stub* stub);
        Matrix<double> H(const Stub* stub);
        Matrix<double> F(const Stub* stub);
        Matrix<double> PxxModel();
        Matrix<double> PddMeas(const Stub* stub);
        std::map<std::string, double> convertParams(std::vector<double> x);
 
    private:
        std::vector<double> mapToVec(std::map<std::string, double> x);
        std::map<std::string, double> vecToMap(std::vector<double> x);
};
#endif
 

