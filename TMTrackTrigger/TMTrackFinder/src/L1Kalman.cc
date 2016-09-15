///=== This is the base class for the Kalman Filter track fit algorithm.

///=== Written by: Sioni Summers

#include "TMTrackTrigger/TMTrackFinder/interface/L1Kalman.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Matrix.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
 
#include <algorithm>
#include <functional>
 
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b){
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::plus<T>());
    return result;
}
template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b){
    assert(a.size() == b.size());
    std::vector<T> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::minus<T>());
    return result;
}
 
L1Kalman::L1Kalman(const Settings* settings, const uint nPar) : TrackFitGeneric(settings){
    nPar_ = nPar;
}
 
L1fittedTrack L1Kalman::fit(const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg){
    const TP* tpa = l1track3D.getMatchedTP();
    if(tpa!=nullptr && getSettings()->debug()==6){
        std::cout << "TP = " << getSettings()->invPtToInvR()*tpa->qOverPt()/2 << "," << tpa->phi0() << ","<<tpa->z0() << ","<< tpa->tanLambda() << std::endl;
    }
    std::vector<const Stub*> stubs = l1track3D.getStubs();
 
    std::vector<double> x0 = seedx(l1track3D);
    Matrix<double> pxxa = seedP(l1track3D);
 
    std::vector<double> xf = x0;
    std::vector<double> xa = x0;
    if (getSettings()->debug()==6) std::cout << "KF = " << xa[0] <<","<< xa[1] <<","<< xa[2] <<"," <<xa[3]<<std::endl;
 
    for(auto &stub : stubs){
        /* Initialise the covariance to a large value to eliminate bias from seed */
        double covMultiplier = 1;
       
        /* Forecast */
        Matrix<double> f = F(stub);
        Matrix<double> h = H(stub);
        xf = f * xa;
        Matrix<double> pxxf = f * pxxa * (f.transpose()) + PxxModel() * covMultiplier;
        Matrix<double> pxdf = pxxf * (h.transpose());
        Matrix<double> pddf = h * pxxf * (h.transpose()) + PddMeas(stub);
        Matrix<double> k = pxdf * (pddf.inverse());
        /* Adjust */
        xa = xf + k * (d(stub) - h * xf);
        std::vector<double> delta = d(stub) - h*xf;
        if (getSettings()->debug()==6) std::cout << "de = " << delta[0] << ","<<delta[1] << std::endl;
        if (getSettings()->debug()==6) std::cout << "KF = " << xa[0] <<","<< xa[1] <<","<< xa[2] <<"," <<xa[3]<<std::endl;
        pxxa = pxxf - k * (pxdf.transpose());
    }
 
    std::map<std::string, double> tp = convertParams(xa);
    return L1fittedTrack(getSettings(), l1track3D, l1track3D.getStubs(), tp["qOverPt"], 0, tp["phi0"], tp["z0"], tp["t"], 0, nPar_, iPhiSec, iEtaReg);
}
 

