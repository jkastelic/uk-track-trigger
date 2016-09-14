#ifndef __KF4ParamsCombV2_H__
#define __KF4ParamsCombV2_H__

#include "TMTrackTrigger/TMTrackFinder/interface/KF4ParamsComb.h"
#include <TMatrixD.h>

class KF4ParamsCombV2 : public KF4ParamsComb{

    public:
	enum PAR_IDS { BETA, Z0P, R0P, RHO0 };  
	enum MEAS_IDS { Z, R };
    public:
	KF4ParamsCombV2(const Settings* settings, const std::string &fitterName );
	~KF4ParamsCombV2(){}
	std::string getParams();

    protected:
	std::map<std::string, double> getTrackParams( const kalmanState *state )const;
	std::vector<double> seedx(const L1track3D& l1track3D)const;
	TMatrixD seedP(const L1track3D& l1track3D)const;
	std::vector<double> d(const Stub* stub )const;
	TMatrixD H(const Stub* stub)const;
	TMatrixD dH(const Stub* stub, const kalmanState *state )const;
	TMatrixD PxxModel( const kalmanState* state, const Stub* stub, unsigned stub_itr )const;
	std::vector<double> ErrMeas(const Stub* stub, std::vector<double> x )const;
	TMatrixD PddMeas(const Stub* stub, const kalmanState *state )const;
	std::map<std::string, double> convertParams(std::vector<double> x)const;
	bool stubBelongs(const Stub* stub, kalmanState& state, std::vector<double> resid)const;
	bool isGoodState( const kalmanState &state )const;

	std::vector<double> residual(const Stub* stub, std::vector<double> &x )const;

};

#endif

