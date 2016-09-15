///=== This is the base class for the Kalman Combinatorial Filter track fit algorithm.
 
///=== Written by: Sioni Summers
 
#ifndef __L1_KALMAN_COMB__
#define __L1_KALMAN_COMB__
 
//#include "TMTrackTrigger/TMTrackFinder/interface/Matrix.h"
#include <TMatrixD.h>
#include "TMTrackTrigger/TMTrackFinder/interface/TrackFitGeneric.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1fittedTrack.h"
#include "TMTrackTrigger/TMTrackFinder/interface/kalmanState.h"
#include <map>
#include <vector>
// #include <fstream>
#include <TString.h>

 
class TH1F;
class TH2F;
class TP; 
class kalmanState;
class L1KalmanComb : public TrackFitGeneric{
 
    public:
        L1KalmanComb(const Settings* settings, const uint nPar, const std::string &fitterName="", const uint nMeas=2 );
 
        virtual ~L1KalmanComb(){}
 
        L1fittedTrack fit(const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg);
	void bookHists();
    protected:
	static  std::map<std::string, double> getTrackParams( const L1KalmanComb *p, const kalmanState *state );
	virtual std::map<std::string, double> getTrackParams( const kalmanState *state )const=0;

	double sectorPhi()const
	{
	    return 2.*M_PI * (0.5 + float(iCurrentPhiSec_)) / float(getSettings()->numPhiSectors()) - M_PI; // Centre of sector in phi
	}
	bool kalmanUpdate( const Stub *stub, kalmanState &state, kalmanState &new_state, const TP *tpa );
	const kalmanState *kalmanUpdate( unsigned nItr, const Stub* stub, const kalmanState &state, const TP *);
	void resetStates();
	const kalmanState *mkState( unsigned nIterations, unsigned layerId, double r, const kalmanState *last_state, 
		const std::vector<double> &x, const TMatrixD &pxx, const Stub* stub, double chi2 );
	//	kalmanState smooth(kalmanState &state);

	virtual std::string getParams()=0;

    protected:
	/* Methods */
	std::vector<double> Hx( const TMatrixD &pH, const std::vector<double> &x )const;
	std::vector<double> Fx( const TMatrixD &pF, const std::vector<double> &x )const;
	TMatrixD HxxH( const TMatrixD &pH, const TMatrixD &xx )const;
	double Chi2( const TMatrixD &dcov, const std::vector<double> &delta, bool debug = false )const;
	TMatrixD GetKalmanMatrix( const TMatrixD &h, const TMatrixD &pxcov, const TMatrixD &dcov )const;
	void GetAdjustedState( const TMatrixD &K, const TMatrixD &pxcov, 
		const std::vector<double> &x, const Stub *stub,  
		std::vector<double> &new_x, TMatrixD &new_xcov )const;


	virtual std::vector<double> seedx(const L1track3D& l1track3D)const=0;
	virtual TMatrixD seedP(const L1track3D& l1track3D)const=0;
	virtual void barrelToEndcap( std::vector<double> &x, TMatrixD &cov_x )const{}
	virtual std::vector<double> d(const Stub* stub )const=0;
	virtual TMatrixD H(const Stub* stub)const=0;
	virtual TMatrixD F(const Stub* stub=0, const kalmanState *state=0 )const=0;
	virtual TMatrixD PxxModel( const kalmanState *state, const Stub *stub, unsigned stub_itr )const=0; 
	virtual std::vector<double> ErrMeas(const Stub* stub, std::vector<double> x )const=0;
	virtual TMatrixD PddMeas(const Stub* stub, const kalmanState *state )const=0;
	virtual bool stubBelongs(const Stub* stub, kalmanState &state, unsigned itr )const=0;

	virtual std::vector<double> residual(const Stub* stub, const std::vector<double> &x )const;
	virtual const kalmanState *updateSeedWithStub( const kalmanState &state, const Stub *stub ){ return 0; }
	virtual bool isGoodState( const kalmanState &state )const{ return true; }

	bool validationGate( const Stub *stub, unsigned stub_itr, const kalmanState &state, double &e2, bool debug = false )const; 
	double calcChi2( unsigned itr, const kalmanState &state )const;
	void printTP( std::ostream &os, const TP *tp )const;


	unsigned getNextLayer( unsigned state_layer, unsigned next_stub_layer );
	std::vector<const Stub *> getNextLayerStubs( const kalmanState *state, std::vector<const Stub *> &stubs, unsigned &next_layer );
	virtual double getRofState( unsigned layerId, const std::vector<double> &xa ) const { return 0;}
	std::vector<const kalmanState *> doKF( unsigned nItr, const std::vector<const kalmanState *> &states, std::vector<const Stub *> stubs, const TP *tpa );

	void fillCandHists( const kalmanState &state, const TP *tpa=0 );
	void fillTrackHists( const kalmanState *state, const TP *tpa, std::vector<const Stub *> &stubs );
	void fillEachNumOfVirtualStubStateHists( unsigned nItr, unsigned nvs0, unsigned nvs1, unsigned nvs2 );
	void fillStepHists( const TP *tpa, unsigned nItr, 
		const TMatrixD &pxxf, const TMatrixD &pxxm, const TMatrixD &pddf,
		const TMatrixD &pddm, const TMatrixD &k, const kalmanState *new_state );

    private:
	unsigned layerContinuity(const Stub* stub, unsigned prevLayerId);
    protected:
	unsigned nPar_;
	unsigned nMeas_;
	std::vector<kalmanState *> state_list_;
	unsigned nIterations_;
	std::vector<double> hkfxmin;
	std::vector<double> hkfxmax;
	std::map<TString, TH1F*> hkfxMap;
	std::vector<double> hxmin;
	std::vector<double> hxmax;
	std::map<TString, TH1F*> hxtMap;
	std::map<TString, TH1F*> hx0Map;
	std::map<TString, TH1F*> hxfMap;
	std::vector<double> hdxmin;
	std::vector<double> hdxmax;
	std::map<TString, TH1F*> hpxxfMap;
	std::map<TString, TH1F*> hpxxaMap;
	std::vector<double> hddmin;
	std::vector<double> hddmax;
	std::map<TString, TH1F*> hpddfMap;
	std::vector<double> hdxModelmin;
	std::vector<double> hdxModelmax;
	std::map<TString, TH1F*> hPxxModelMap;
	std::vector<double> hddMeasmin;
	std::vector<double> hddMeasmax;
	std::map<TString, TH1F*> hPddMeasMap;
	std::vector<double> hresmin;
	std::vector<double> hresmax;
	std::map<TString, TH1F*> hresMap;
	std::map<TString, TH1F*> hkMap;
	std::vector<double> hxaxtmin;
	std::vector<double> hxaxtmax;
	std::map<TString, TH1F*> hxaxtMap;
	std::map<TString, TH1F*> hNumStatesItrMap;
	TH1F*           hnmergeStub_;

	TH2F* hTrackStubEtaLayer;
	double hchi2min;
	double hchi2max;
	std::map<TString, TH1F*> hchi2Map;

	unsigned maxNfitForDump_;
	bool     dump_;
	unsigned int      iCurrentPhiSec_;
	unsigned int      iCurrentEtaReg_;
};
#endif




