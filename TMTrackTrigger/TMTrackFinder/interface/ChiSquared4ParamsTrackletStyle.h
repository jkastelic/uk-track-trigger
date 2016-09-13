#ifndef __CHI_SQUARED_4_PARAMS_TRACKLET_STYLE__
#define __CHI_SQUARED_4_PARAMS_TRACKLET_STYLE__
 
#include "TMTrackTrigger/TMTrackFinder/interface/L1ChiSquared.h"
 
class ChiSquared4ParamsTrackletStyle : public L1ChiSquared{
 
public:
    ChiSquared4ParamsTrackletStyle(const Settings* settings, const uint nPar);
 
    ~ChiSquared4ParamsTrackletStyle(){}
 
    std::string getParams();
 
protected:
    std::vector<double> seed(const L1track3D& l1track3D);
    std::vector<double> residuals(std::vector<double> x);
    Matrix<double> D(std::vector<double> x);
    Matrix<double> Vinv();
    std::map<std::string, double> convertParams(std::vector<double> x);
 
private:
    std::vector<double> mapToVec(std::map<std::string, double> x);
    std::map<std::string, double> vecToMap(std::vector<double> x);
};
 
#endif

