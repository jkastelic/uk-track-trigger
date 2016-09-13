#ifndef __KALMAN_STATE__
#define __KALMAN_STATE__
 
//#include "TMTrackTrigger/TMTrackFinder/interface/Matrix.h"
#include <TMatrixD.h>
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1KalmanComb.h"
#include <map>

class L1KalmanComb;
class kalmanState;
typedef std::map<std::string, double> (*GET_TRACK_PARAMS)( const L1KalmanComb *p, const kalmanState *state );
 
class kalmanState{
    public:
	kalmanState();
	kalmanState( unsigned nIterations, unsigned layerId, const kalmanState *last_state, const std::vector<double> &x, const TMatrixD &pxx, const Stub* stub, double chi2, 
		L1KalmanComb *fitter, GET_TRACK_PARAMS f );
	kalmanState(const kalmanState &p);
	~kalmanState(){}

	kalmanState & operator=( const kalmanState &other );

	unsigned         nIterations()const{ return nIterations_; }
	unsigned             layerId()const{ return     layerId_; }
	bool                  barrel()const{ return      barrel_; }
	double                     r()const{ return           r_; }
	double                     z()const{ return           z_; }
	const kalmanState *last_state()const{ return  last_state_; }
	std::vector<double>       xa()const{ return          xa_; }
	TMatrixD                pxxa()const{ return        pxxa_; }
	const Stub*             stub()const{ return        stub_; }
	double                  chi2()const{ return        chi2_; }
	unsigned              nStubs()const{ return      n_stubs_; }
	unsigned       nVirtualStubs()const{ return n_virtual_stubs_; }
	unsigned         nStubLayers()const{ return n_stub_layers_; }
	bool                    good( const TP *tp )const;
	double           reducedChi2()const;
	const kalmanState *last_update_state()const;
	std::vector<const Stub *>      stubs()const;
	L1KalmanComb      *fitter()const{ return fitter_; }
	GET_TRACK_PARAMS fXtoTrackParams()const{ return fXtoTrackParams_; };

	static bool orderReducedChi2(const kalmanState *left, const kalmanState *right);
	static bool order(const kalmanState *left, const kalmanState *right);
	void dump( ostream &os, const TP *tp=0, bool all=0 )const;
	void setChi2( double p ){ chi2_ = p; }

    private:
	unsigned             nIterations_;
	unsigned                 layerId_;
	double                         r_;
	const kalmanState    *last_state_;
	std::vector<double>           xa_;
	TMatrixD                    pxxa_;
	const Stub                 *stub_;
	double                      chi2_;
	unsigned                 n_stubs_;
	unsigned         n_virtual_stubs_;
	unsigned           n_stub_layers_;
	L1KalmanComb             *fitter_;
	GET_TRACK_PARAMS fXtoTrackParams_;
	bool                      barrel_;
	double                         z_;

};
#endif


