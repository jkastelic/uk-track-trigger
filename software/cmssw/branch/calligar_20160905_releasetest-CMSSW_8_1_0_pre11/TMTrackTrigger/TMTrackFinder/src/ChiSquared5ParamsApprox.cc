#include "TMTrackTrigger/TMTrackFinder/interface/ChiSquared5ParamsApprox.h"
 
ChiSquared5ParamsApprox::ChiSquared5ParamsApprox(const Settings* settings, const uint nPar) : L1ChiSquared(settings, nPar){
    //parameterStream_ << "4Params_TrackletStyle_MCTruthSeed";
    //configParameters_ = (lsStr_.str());
 
}
 
std::map<std::string, double> ChiSquared5ParamsApprox::vecToMap(std::vector<double> x){
    // Convert a vector of track parameters to a labelled map for ease of use
    std::map<std::string, double> result;
    result["2r1Inv"] = x[0];
    result["phi0"] = x[1];
    result["d1"] = x[2];
    result["z0"] = x[3];
    result["t"] = x[4];
    return result;
}
 
std::vector<double> ChiSquared5ParamsApprox::mapToVec(std::map<std::string, double> x){
    // Conevrt the map of labelled track parameters to a vector (in correct order)
    std::vector<double> result;
    result.resize(5);
    result[0] = x["2r1Inv"];
    result[1] = x["phi0"];
    result[2] = x["d1"];
    result[3] = x["z0"];
    result[4] = x["t"];
    return result;
}
 
std::vector<double> ChiSquared5ParamsApprox::seed(const L1track3D& l1track3D){
    /* Cheat by using MC trutth to initialize helix parameters. Useful to check if conevrgence is the problem */
    std::map<std::string, double> x;
    double d0 = l1track3D.d0();
    double r0 = 1/(getSettings()->invPtToInvR() * l1track3D.qOverPt());
    x["2r1Inv"] = 1/(2*(r0 + d0));
    x["phi0"] = l1track3D.phi0();
    x["d1"] = d0 * (1 - d0/(2*(r0 + d0)));
    x["z0"] = l1track3D.z0();
    x["t"] = l1track3D.tanLambda();
    return mapToVec(x);
}
 
Matrix<double> ChiSquared5ParamsApprox::D(std::vector<double> x){
    Matrix<double> D(2 * stubs_.size(), nPar_, 0.0); // Empty matrix
    int j = 0;
    std::map<std::string, double> y = vecToMap(x); // Get the track params by label
    double rInv = y["2r1Inv"];
    double phi0 = y["phi0"];
    double d0 = y["d1"];
    double t = y["t"];
    double z0 = y["z0"];
    for(unsigned i = 0; i < stubs_.size(); i++){
        double ri=stubs_[i]->r();
        double zi=stubs_[i]->z();
   
        D(j, 0) = ri;
        D(j, 1) = 1;
        D(j, 2) = 1/ri;
        //D(j, 3) = 0;
 
        j++;
        //D(j, 0)
        //D(j, 1)
        D(j, 3) = 1;
        D(j, 4) = ri;
        j++;
    }
    return D;
}
 
Matrix<double> ChiSquared5ParamsApprox::Vinv(){
    Matrix<double> Vinv(2*stubs_.size(), 2*stubs_.size(), 0.0);
    for(unsigned i = 0; i < stubs_.size(); i++){
        Vinv(2*i, 2*i) = 1/stubs_[i]->sigmaX();
        Vinv(2*i + 1, 2*i + 1) = 1/stubs_[i]->sigmaZ();
 
    }
    return Vinv;
}
 
std::vector<double> ChiSquared5ParamsApprox::residuals(std::vector<double> x) {
 
    unsigned int n=stubs_.size();
   
    std::vector<double> delta;
    delta.resize(2*n);
 
    std::map<std::string, double> trackParams = vecToMap(x); // Get the track params by label
    double rInv = trackParams["2r1Inv"];
    double phi0 = trackParams["phi0"];
    double d1 = trackParams["d1"];
    double t = trackParams["t"];
    double z0 = trackParams["z0"];
 
    double chisq=0.0;
 
    unsigned int j=0;
 
    bool print=false;
 
    if (print) std::cout << "Residuals ("<<chisq<<") ["<<getSettings()->invPtToInvR()/rInv<<"]: ";

    float largestresid=-1.0;
    int ilargestresid=-1;
 
    for(unsigned int i=0;i<n;i++) {
        double ri=stubs_[i]->r();
        double zi=stubs_[i]->z();
        double phii=stubs_[i]->phi();
        const double sigmax=stubs_[i]->sigmaX();
        const double sigmaz=stubs_[i]->sigmaZ();
 
        delta[j++] = phii - (phi0 + ri*rInv + d1/ri);
        delta[j++] = zi - (z0 + ri*t);
       
        if (print) std::cout << delta[j-2]<<" "<<delta[j-1]<<" ";
 
        chisq+=delta[j-2]*delta[j-2]+delta[j-1]*delta[j-1];

        if (fabs(delta[j-2])>largestresid) {
          largestresid=fabs(delta[j-2]);
          ilargestresid=i;
        }
 
        if (fabs(delta[j-1])>largestresid) {
          largestresid=fabs(delta[j-1]);
          ilargestresid=i;
        }
 
    if (print) std::cout <<" ("<<chisq<<")"<< std::endl;
    }

    largestresid_ = largestresid;
    ilargestresid_ = ilargestresid;
    return delta;
}
 
std::map<std::string, double> ChiSquared5ParamsApprox::convertParams(std::vector<double> x){
    std::map<std::string, double> y = vecToMap(x); // Get track parameters by label
    std::map<std::string, double> result;
    result["qOverPt"] = 2 * y["2rInv"] / getSettings()->invPtToInvR();
    result["phi0"] = y["phi0"];
    result["z0"] = y["z0"];
    result["t"] = y["t"];
    return result;
}
 
std::string ChiSquared5ParamsApprox::getParams(){
    return "ChiSquared5ParamsApprox";
}
 

