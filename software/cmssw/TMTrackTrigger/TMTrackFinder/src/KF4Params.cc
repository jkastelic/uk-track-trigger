///=== This is the Kalman Filter for 4 helix parameters track fit algorithm.

///=== Written by: Sioni Summers

#include "TMTrackTrigger/TMTrackFinder/interface/KF4Params.h"
 
KF4Params::KF4Params(const Settings* settings, const uint nPar) : L1Kalman(settings, nPar){
}
 
std::vector<double> KF4Params::mapToVec(std::map<std::string, double> x){
    std::vector<double> y;
    y.resize(4);
    y[0] = x["rInv"];
    y[1] = x["phi0"];
    y[2] = x["z0"];
    y[3] = x["t"];
    return y;
}
 
std::map<std::string, double> KF4Params::vecToMap(std::vector<double> x){
    std::map<std::string, double> y;
    y["rInv"] = x[0];
    y["phi0"] = x[1];
    y["z0"] = x[2];
    y["t"] = x[3];
    return y;
}
 
Matrix<double> KF4Params::H(const Stub* stub){
    Matrix<double> h(2, 4, 0.0);
    h(0,0) = -stub->r();
    h(0,1) = 1;
    h(1,2) = 1;
    h(1,3) = stub->r();
    return h;
}
 
std::vector<double> KF4Params::seedx(const L1track3D& l1track3D){
    std::map<std::string, double> x;
    x["rInv"] = getSettings()->invPtToInvR() * l1track3D.qOverPt()/2;
    x["phi0"] = l1track3D.phi0();
    x["z0"] = l1track3D.z0();
    x["t"] = l1track3D.tanLambda();
    return mapToVec(x);
}
 
Matrix<double> KF4Params::seedP(const L1track3D& l1track3D){
    Matrix<double> p(4,4,0);
    for(int n = 0; n < 4; n++)
        p(n,n) = 1024;
    return p;
}
 
Matrix<double> KF4Params::F(const Stub* stub){
    Matrix<double> F(4,4,0.0);
    for(int n = 0; n < 4; n++)
        F(n, n) = 1;
    return F;
}
 
std::vector<double> KF4Params::d(const Stub* stub){
    std::vector<double> meas;
    meas.resize(2);
    meas[0] = stub->phi();
    meas[1] = stub->z();
    return meas;
}
 
Matrix<double> KF4Params::PddMeas(const Stub* stub){
    Matrix<double> p(2,2,0);
    p(0,0) = stub->sigmaX();
    p(1,1) = stub->sigmaZ();
    return p;
}
 
Matrix<double> KF4Params::PxxModel(){
    // TODO DO THIS PROPERLY
    Matrix<double> p(4,4,0);
    for(int n = 0; n < 4; n++)
        p(n,n) = 0.1;
    return p;
}
 
std::map<std::string, double> KF4Params::convertParams(std::vector<double> x){
    std::map<std::string, double> y = vecToMap(x);
    std::map<std::string, double> z;
    z["qOverPt"] = (2*y["rInv"]) / getSettings()->invPtToInvR();
    z["phi0"] = y["phi0"];
    z["z0"] = y["z0"];
    z["t"] = y["t"];
    return z;
}
 
std::string KF4Params::getParams(){
    return "KF4Params";
}
 

