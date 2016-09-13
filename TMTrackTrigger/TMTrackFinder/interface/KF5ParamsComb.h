#ifndef __KF5ParamsComb_H__
#define __KF5ParamsComb_H__

#include "TMTrackTrigger/TMTrackFinder/interface/L1KalmanComb.h"
#include <TMatrixD.h>

class KF5ParamsComb : public L1KalmanComb{

    public:
	enum PAR_IDS { INV2R, PHI0, Z0, T, D0 }; 
	enum MEAS_IDS { PHI, Z };

    public:
	KF5ParamsComb(const Settings* settings, const string &fitterName );
	~KF5ParamsComb(){}
	std::string getParams();

    protected:

	std::map<std::string, double> getTrackParams(const kalmanState *state )const;
	std::vector<double> seedx(const L1track3D& l1track3D)const;
	TMatrixD seedP(const L1track3D& l1track3D)const;
	std::vector<double> d(const Stub* stub )const;
	TMatrixD H(const Stub* stub)const;
	TMatrixD F(const Stub* stub = 0, const kalmanState *state = 0 )const;
	TMatrixD PxxModel( const kalmanState *state, const Stub *stub, unsigned stub_itr )const; 
	std::vector<double> ErrMeas(const Stub* stub, std::vector<double> x )const;
	TMatrixD PddMeas(const Stub* stub, const kalmanState *state )const;
	std::map<std::string, double> convertParams(std::vector<double> x)const;
	bool stubBelongs(const Stub* stub, kalmanState& state, unsigned itr )const;

	std::vector<double> residual(const Stub* stub, const std::vector<double> &x )const;
	const kalmanState *updateSeedWithStub( const kalmanState &state, const Stub *stub );
	bool isGoodState( const kalmanState &state )const;

	double getRofState( unsigned layerId, const vector<double> &xa )const;
	TMatrixD dH(const Stub* stub)const;

};

#endif

