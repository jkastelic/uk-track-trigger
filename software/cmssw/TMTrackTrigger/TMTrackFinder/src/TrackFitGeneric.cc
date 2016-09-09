///=== This is the base class for all the track fit algorithms

///=== Written by: Alexander D. Morton and Sioni Summers

#include "TMTrackTrigger/TMTrackFinder/interface/TrackFitGeneric.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Utility.h"
#include "TMTrackTrigger/TMTrackFinder/interface/ChiSquared4ParamsTrackletStyle.h"
#include "TMTrackTrigger/TMTrackFinder/interface/ChiSquared4ParamsApprox.h"
#include "TMTrackTrigger/TMTrackFinder/interface/ChiSquared5ParamsApprox.h"
#include "TMTrackTrigger/TMTrackFinder/interface/KF4Params.h"
#include "TMTrackTrigger/TMTrackFinder/interface/KF4ParamsComb.h"
#include "TMTrackTrigger/TMTrackFinder/interface/KF4ParamsCombV2.h"
//#include "TMTrackTrigger/TMTrackFinder/interface/KF4ParamsCombV3.h"
#include "TMTrackTrigger/TMTrackFinder/interface/KF5ParamsComb.h"
//#include "TMTrackTrigger/TMTrackFinder/interface/KF4ParamsCombIV.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TrackFitLinearAlgo.h"
#include "TMTrackTrigger/TMTrackFinder/interface/LinearRegression.h"
#include <map> 
#include <new>
 
//=== Set configuration parameters.
 
TrackFitGeneric::TrackFitGeneric( const Settings* settings, const string &fitterName ) : settings_(settings), fitterName_(fitterName), nDupStubs_(0) {
}
 
 
//=== Fit a track candidate obtained from the Hough Transform.
//=== Specify which phi sector and eta region it is in.
 
L1fittedTrack TrackFitGeneric::fit(const L1track3D& l1track3D,  unsigned int iPhiSec, unsigned int iEtaReg) {
  return L1fittedTrack (settings_, l1track3D, l1track3D.getStubs(), 0, 0, 0, 0, 0, 999999., 0, iPhiSec, iEtaReg);
}
 
/*std::auto_ptr<TrackFitGeneric> TrackFitGeneric::create(std::string fitter, const Settings* settings){
    if(fitter.compare("ChiSquared4ParamsTrackletStyle") == 0){
        return std::auto_ptr<TrackFitGeneric>(new ChiSquared4ParamsTrackletStyle(settings, 4));
    }else if(fitter.compare("TrackFitLinearAlgo4") == 0){
        return std::auto_ptr<TrackFitGeneric>(new TrackFitLinearAlgo(settings, 4));
    }else if(fitter.compare("TrackFitLinearAlgo5") == 0){
        return std::auto_ptr<TrackFitGeneric>(new TrackFitLinearAlgo(settings, 5));
    }else if(fitter.compare("short")==0){
        return std::auto_ptr<TrackFitGeneric>(new ChiSquared4ParamsTrackletStyle(settings, 4));
    }else if(fitter.compare("ChiSquared4ParamsApprox")==0){
        return std::auto_ptr<TrackFitGeneric>(new ChiSquared4ParamsApprox(settings, 4));
    }else if(fitter.compare("ChiSquared5ParamsApprox")==0){
        return std::auto_ptr<TrackFitGeneric>(new ChiSquared5ParamsApprox(settings, 5));
    }else if(fitter.compare("KF4Params")==0){
        return std::auto_ptr<TrackFitGeneric>(new KF4Params(settings, 4));
    }else if(fitter.compare("KF4ParamsComb")==0){
        return std::auto_ptr<TrackFitGeneric>(new KF4ParamsComb(settings, 4));
    }else if(fitter.compare("LinearRegression")==0){
        return std::auto_ptr<TrackFitGeneric>(new LinearRegression(settings, 4));
    }else{
        std::cout << "Requested unknown fitter " << fitter << ", using TrackFitLinearAlgo instead" << std::endl;
        return std::auto_ptr<TrackFitGeneric>(new TrackFitLinearAlgo(settings, 4));
    }
}
*/ 

TrackFitGeneric* TrackFitGeneric::create(std::string fitter, const Settings* settings){
    if(fitter.compare("ChiSquared4ParamsTrackletStyle") == 0){
	return new ChiSquared4ParamsTrackletStyle(settings, 4);
    }else if(fitter.compare("TrackFitLinearAlgo4") == 0){
	return new TrackFitLinearAlgo(settings, 4);
    }else if(fitter.compare("TrackFitLinearAlgo5") == 0){
	return new TrackFitLinearAlgo(settings, 5);
    }else if(fitter.compare("short")==0){
	return new ChiSquared4ParamsTrackletStyle(settings, 4);
    }else if(fitter.compare("ChiSquared4ParamsApprox")==0){
	return new ChiSquared4ParamsApprox(settings, 4);
    }else if(fitter.compare("ChiSquared5ParamsApprox")==0){
	return new ChiSquared5ParamsApprox(settings, 5);
    }else if(fitter.compare("KF4Params")==0){
	return new KF4Params(settings, 4);
    }else if(fitter.compare("KF4ParamsComb")==0){
	return new KF4ParamsComb(settings, 4, fitter );
    }else if(fitter.compare("KF4ParamsCombV2")==0){
	return new KF4ParamsCombV2(settings, fitter );
	//    }else if(fitter.compare("KF4ParamsCombV3")==0){
	//	return new KF4ParamsCombV3(settings, fitter );
    }else if(fitter.compare("KF5ParamsComb")==0){
	return new KF5ParamsComb(settings, fitter );
	//    }else if(fitter.compare("KF4ParamsCombIV")==0){
	//	return new KF4ParamsCombIV(settings, fitter );
    }else if(fitter.compare("LinearRegression")==0){
        return new LinearRegression(settings, 4);
    }else{
	std::cout << "Requested unknown fitter " << fitter << ", using TrackFitLinearAlgo instead" << std::endl;
	return new TrackFitLinearAlgo(settings, 4);
    }
} 

