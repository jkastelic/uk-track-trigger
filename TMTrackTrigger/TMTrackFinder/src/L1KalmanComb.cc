///=== This is the base class for the Kalman Combinatorial Filter track fit algorithm.

///=== Written by: Sioni Summers

#include "TMTrackTrigger/TMTrackFinder/interface/L1KalmanComb.h"
//#include "TMTrackTrigger/TMTrackFinder/interface/Matrix.h"
#include <TMatrixD.h> 
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
#include "TMTrackTrigger/TMTrackFinder/interface/kalmanState.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <algorithm>
#include <functional>
#include <fstream>
#include <iomanip>
#include <TH2F.h>
//#define CKF_DEBUG

using namespace std;

static double wrapRadian( double t ){

    if( t > 0 ){
	while( t > M_PI ) t-= 2*M_PI; 
    }
    else{
	while( t < - M_PI ) t+= 2*M_PI; 
    }
    return t;
}

class StubWithValidation{

    public :
	StubWithValidation(): stub_(0), e2_(0){}
	StubWithValidation( const Stub *stub, double e2 ):stub_(stub), e2_(e2){}
	~StubWithValidation(){}
	const Stub *stub()const{ return stub_; }
	double e2()const{ return e2_; }
	void set_e2( double e2 ){ e2_ = e2; }
    private:
	const Stub *stub_;
	double e2_;

};
/*
static bool order_with_validation( const StubWithValidation &a, const StubWithValidation &b ){
    return ( a.e2() < b.e2() );
}
*/

unsigned LayerId[16] = { 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 21, 22, 23, 24, 25 };

static bool orderStubsByLayer(const Stub* a, const Stub* b){
    return (a->layerId() < b->layerId());
}

static bool orderStubsByZ(const Stub* a, const Stub* b){
    return (a->z() < b->z());
}
static bool orderStubsByR(const Stub* a, const Stub* b){
    return (a->r() < b->r());
}
static bool isOverlap( const Stub* a, const Stub*b ){
    if( a->layerId() != b->layerId() ) return false;
    if( a->layerId() < 7 ){
	if( fabs( b->z() - a->z() ) > a->stripLength() || fabs( b->phi() - a->phi() ) > a->stripPitch() / a->r()  ) return false;
    }
    else{
	if( fabs( b->r() - a->r() ) > a->stripLength() || fabs( b->phi() - a->phi() ) > a->stripPitch() / a->r() ) return false;
    }
    return true;
}


std::map<std::string, double> L1KalmanComb::getTrackParams( const L1KalmanComb *p, const kalmanState *state )
{
    return p->getTrackParams( state );
}
std::vector<double> L1KalmanComb::Hx( const TMatrixD &pH, const std::vector<double> &x )const
{
    std::vector<double> m( (unsigned) pH.GetNrows() );
    if( pH.GetNcols() != (int) x.size() ) { std::cerr << "Hx() : H and x have different dimensions" << std::endl; }
    else{

	for( int i=0; i < pH.GetNcols(); i++ ){ 
	    for( int j=0; j < pH.GetNrows(); j++ ){ 
		m.at(j) += pH(j,i) * x.at(i);
	    }
	}
    }
    return m;
}
std::vector<double> L1KalmanComb::Fx( const TMatrixD &pF, const std::vector<double> &x )const
{
    return Hx( pF, x );
}
TMatrixD L1KalmanComb::HxxH( const TMatrixD &pH, const TMatrixD &xx )const
{
    int nd = (unsigned) pH.GetNrows(); 
    TMatrixD tmp(nd,nPar_);
    TMatrixD mHxxH(nd,nd);
    if( pH.GetNcols() != xx.GetNcols() || pH.GetNcols() != xx.GetNrows() ) { std::cerr << "HxxH() : H and xx have different dimensions" << std::endl; }
    else{

	for( int i=0; i < pH.GetNrows(); i++ ){ 
	    for( int j=0; j < xx.GetNrows(); j++ ){ 
		for( int k=0; k < xx.GetNcols(); k++ ){ 
		    tmp(i,k) += pH(i,j) * xx(j,k);
		}
	    }
	}
	for( int i=0; i < tmp.GetNrows(); i++ ){ 
	    for( int j=0; j < pH.GetNcols(); j++ ){ 
		for( int k=0; k < pH.GetNrows(); k++ ){ 
		    mHxxH(i,k) += tmp(i,j) * pH(k,j); 
		}
	    }
	}
    }
    return mHxxH;

}
double L1KalmanComb::Chi2( const TMatrixD &dcov, const std::vector<double> &delta, bool debug )const
{
      if( dcov.Determinant() == 0 ) return 999;

    TMatrixD dcovi( dcov );
    dcovi.Invert();

    vector<double> tmp(2,0);
    for( int i=0; i < dcovi.GetNrows(); i++ ){ 
	for( int j=0; j < dcovi.GetNcols(); j++ ){ 
	    tmp.at(j) += delta.at(i) * dcovi(i,j); 
	}
    }
    double chi2(0);
    for( int j=0; j < 2; j++ ){ 
	chi2 += tmp.at(j) * delta.at(j);
    }

    if( debug ){
	cout << "CHI SQUARE OUTPUT" << endl;
	cout << "cov" << endl;
	dcov.Print();
	cout << "cov inv" << endl;
	dcovi.Print();
	for( unsigned i=0; i < delta.size(); i++ ) cout << delta.at(i) << " ";
	cout << endl;
    }
    return chi2;
}
TMatrixD L1KalmanComb::GetKalmanMatrix( const TMatrixD &h, const TMatrixD &pxcov, const TMatrixD &dcov )const
{
  
    TMatrixD pxcovht(pxcov.GetNrows(),2);
    for( int i=0; i<pxcov.GetNrows(); i++ ){
	for( int j=0; j<pxcov.GetNcols(); j++ ){
	    for( int k=0; k<h.GetNrows(); k++ ){
		pxcovht(i,k) += pxcov(i,j) * h(k,j);
	    }
	}
    }
    if( getSettings()->kalmanDebugLevel() >= 4 ){
	cout << "pxcovht" << endl;
	pxcovht.Print();
    }

    TMatrixD tmp(dcov.GetNrows(), dcov.GetNcols() );
    TMatrixD hxxh = HxxH( h, pxcov );
    tmp = dcov + hxxh; 

    if( getSettings()->kalmanDebugLevel() >= 4 ){
	cout << "hxxh" << endl;
	hxxh.Print();
	cout << "dcov + hxxh " << endl;
	tmp.Print();
    }

    TMatrixD K( pxcovht.GetNrows(), tmp.GetNcols() );

    if(tmp.Determinant() == 0 ) return K; 
    tmp.Invert();

    for( int i=0; i<pxcovht.GetNrows(); i++ ){
	for( int j=0; j<pxcovht.GetNcols(); j++ ){
	    for( int k=0; k<tmp.GetNcols(); k++ ){
		K(i,k)+=pxcovht(i,j)*tmp(j,k);
	    }
	}
    }
    return K;
}

void L1KalmanComb::GetAdjustedState( const TMatrixD &K, const TMatrixD &pxcov, 
	const std::vector<double> &x, const Stub *stub, 
	std::vector<double> &new_x, TMatrixD &new_xcov )const
{
    TMatrixD h = H(stub);
    std::vector<double> m = d(stub );

    for( unsigned i=0; i < new_x.size(); i++ ){
	new_x.at(i) = 0;
    }

    std::vector<double> tmpv(m.size(), 0 );
    for( int i=0; i < h.GetNrows(); i++ ){
	tmpv.at(i) += m.at(i); 
	for( int j=0; j < h.GetNcols(); j++ ){
	    tmpv.at(i) += -1. * h(i,j) * x.at(j); 
	}
    }
    for( int i=0; i < K.GetNrows(); i++ ){
	new_x.at(i) += x.at(i);
	for( int j=0; j < K.GetNcols(); j++ ){
	    new_x.at(i) += K(i,j) * tmpv.at(j);
	}
    }

    TMatrixD tmp(K.GetNrows(), h.GetNcols() );
    for( int i=0; i< K.GetNrows(); i++ ){
	tmp(i,i) = 1;
    }
    for( int i=0; i< K.GetNrows(); i++ ){
	for( int j=0; j< K.GetNcols(); j++ ){
	    for( int k=0; k< h.GetNcols(); k++ ){
		tmp(i,k) += -1 * K(i,j) * h(j,k);
	    }
	}
    }
    new_xcov.Clear();
    new_xcov.ResizeTo(pxcov.GetNrows(), pxcov.GetNcols());
    for( int i=0; i< tmp.GetNrows(); i++ ){
	for( int j=0; j< tmp.GetNcols(); j++ ){
	    for( int k=0; k< pxcov.GetNcols(); k++ ){
		new_xcov(i,k) += tmp(i,j) * pxcov(j,k);
	    }
	}
    }
}

void L1KalmanComb::printTP( std::ostream &os, const TP *tp )const{
  
    std::map<std::string, double> tp_x;
    bool useForAlgEff(false);
    if( tp ){
	useForAlgEff = tp->useForAlgEff();
	tp_x["qOverPt"] = tp->qOverPt();
	tp_x["phi0"] = tp->phi0();
	tp_x["z0"] = tp->z0();
	tp_x["t"] = tp->tanLambda();
	tp_x["d0"] = tp->d0();
    }
    if( tp ){
	os << "\tTP index = " << tp->index() << " useForAlgEff = " << useForAlgEff << " ";
	os << "\tpT, eta = " << tp->pt() << ", " << tp->eta() << " ";
	os << "\tqOver2R0 = " << tp->qOverPt() * getSettings()->invPtToInvR() * 0.5 << " "; 
	for( auto pair : tp_x ){
	    os << pair.first << ":" << pair.second << ", "; 
	}
    }
    else{
	os << "\tTP index = "; 
    }
    os << endl;
}
static void printStubLayers( std::ostream &os, std::vector<const Stub *> &stubs ){

    if( stubs.size() == 0 ) os << "stub layers = []" << std::endl;
    else{
	os << "stub layers = [ ";
	for( unsigned i=0; i<stubs.size()-1; i++ ) os << stubs[i]->layerId() << ", ";
	os << stubs.back()->layerId() << " ]" << endl;
    }
}
static void printStub( std::ostream &os, const Stub * stub ){
      os << "stub ";
    os << "[" << stub << "] "; 
    os << "layerId : " << stub->layerId() << " ";
    os << "index : " << stub->index() << " ";
    os  << "[r,phi,z] = ";
    os << "[" << stub->r() << ", " << stub->phi() << ", " << stub->z() << "] ";
    os << " assoc TP indices = [ "; 
    std::set<const TP*> tps = stub->assocTPs();
    for( auto tp : tps ) os << tp->index() << " "; 
    os << "] ";
    os << endl;

}
static void printStubs( std::ostream &os, std::vector<const Stub *> &stubs ){

    for( auto &stub : stubs ){
	printStub( os, stub );
    }
}
/*
   static void printStubAssociatedTPs( std::ostream &os, std::vector<const Stub *> &stubs ){

   for( unsigned i=0; i<stubs.size(); i++ ){
   std::set<const TP*> tps = stubs[i]->assocTPs();
   os << "stub TP indices = [ ";
   for( auto tp : tps ) os << tp->index() << " "; 
   os << "] ";
   }
   os << endl;
   }
   */

L1KalmanComb::L1KalmanComb(const Settings* settings, const uint nPar, const string &fitterName, const uint nMeas ) : TrackFitGeneric(settings, fitterName ){
      
    nPar_ = nPar;
    nMeas_ = nMeas;
    hkfxmin = vector<double>( nPar_, -1 );
    hkfxmax = vector<double>( nPar_,  1 );
    hxmin = vector<double>( nPar_, -1 );
    hxmax = vector<double>( nPar_,  1 );
    hxmin[0] = -0.05;
    hxmax[0] = +0.05;
    hxmin[1] = -3.2;
    hxmax[1] = +3.2;
    hxmin[2] = -20;
    hxmax[2] = +20;
    hxmin[3] = -10;
    hxmax[3] = +10;
    hxmin[4] = -5;
    hxmax[4] = +5;


    hdxmin = vector<double>( nPar_, -1e-4 );
    hdxmax = vector<double>( nPar_,  1e-4 );

    hdxModelmin = vector<double>( nPar_, -1e-6 );
    hdxModelmax = vector<double>( nPar_,  1e-6 );

    hddmin = vector<double>( nPar_, -1e-2 );
    hddmax = vector<double>( nPar_,  1e-2 );

    hddMeasmin = vector<double>( 2, -1e-3 );
    hddMeasmax = vector<double>( 2,  1e-3 );

    hresmin = vector<double>( 2, -1e-2 );
    hresmax = vector<double>( 2,  1e-2 );

    hxaxtmin = vector<double>( nPar_, -1 );
    hxaxtmax = vector<double>( nPar_,  1 );

    hchi2min = 0; 
    hchi2max = 50; 

    maxNfitForDump_ = 10; 
    dump_ = false; 

}

L1fittedTrack L1KalmanComb::fit(const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg){
    
    iCurrentPhiSec_ = iPhiSec;
    iCurrentEtaReg_ = iEtaReg;
    resetStates();

    //TP
    const TP* tpa(0);
    if( l1track3D.getMatchedTP() ){
	tpa = l1track3D.getMatchedTP();
    }
    /*
       if( tpa == 0 || !tpa->useForAlgEff() ){ 
       L1fittedTrack returnTrk(getSettings(), l1track3D, l1track3D.getStubs(), l1track3D.qOverPt(), 0, l1track3D.phi0(), 0, 0, 9999, nPar_, iPhiSec, iEtaReg, false);
       return returnTrk;
       }
       */

    //dump flag
    static unsigned nthFit(0);
    nthFit++;
    if( getSettings()->kalmanDebugLevel() > 2 && nthFit <= maxNfitForDump_ ){
	if( tpa ) dump_ = true; 
	else dump_ = false;
    }
    else dump_ = false;

    //stub list from L1track3D, sorted in layer order
    std::vector<const Stub*> stubs = l1track3D.getStubs();
    sort(stubs.begin(), stubs.end(), orderStubsByLayer); 

    for(unsigned i=0; i < stubs.size(); i++ ){
	const Stub *stub_a = stubs.at(i);
	for(unsigned j=i+1; j < stubs.size(); j++ ){
	    const Stub *stub_b = stubs.at(j);
	    if( stub_a->r() == stub_b->r() && stub_a->phi() == stub_b->phi() && stub_a->z() == stub_b->z() ){
		stubs.erase( stubs.begin() + j ); 
		nDupStubs_ ++;
		j--;
	    }
	}
    }


    //seed
    std::vector<double> x0 = seedx(l1track3D);
    TMatrixD pxx0 = seedP(l1track3D);
    std::vector<const kalmanState *> states;

    const kalmanState *state0 = mkState( 0, 0, 0, 0, x0, pxx0, 0, 0 );
    states.push_back( state0 );

    //fill histograms for the track informations
    if( getSettings()->kalmanFillInternalHists() ) fillTrackHists( state0, tpa, stubs );

    //track information dump
    if( getSettings()->kalmanDebugLevel() >= 1 ){

	std::cout << "===============================================================================" << endl;
	std::cout << "Track Finding candidate in [phi_sec, eta_reg] = [" << iPhiSec << ", " << iEtaReg << "]" << std::endl;
	printTP( cout, tpa );
	printStubLayers( cout, stubs );
	printStubs( cout, stubs );

	cout << "The seed state : ";
	state0->dump( cout, tpa );
    }

    //Kalman Filter
    std::vector<const kalmanState *> last_states = doKF( 1, states, stubs, tpa );
    //sort the candidate states in # of layer, more stubs come first
    sort( last_states.begin(), last_states.end(), kalmanState::order);

    //final states dump
    if( getSettings()->kalmanDebugLevel() >= 1 ){
	cout << "------------------------------------" << endl;
	cout << "# of final states : " << last_states.size() << endl;
	if( getSettings()->kalmanDebugLevel() >= 2 ){
	    if( last_states.size() > 50 ) {
		printStubs( cout, stubs );
		//		printStubAssociatedTPs( cout, stubs );
	    }
	    for( auto &last_state : last_states ) last_state->dump( cout, tpa );
	}
	cout << "------------------------------------" << endl;
    }

    //single state selection, update the candidate with smaller chisquare for the state with enough stubs.
    double cand_chi2(0);
    unsigned last_nstub_layers(0);
    const kalmanState *cand(0);
    std::vector<const kalmanState *> best_cands;
    std::vector<const kalmanState *>::iterator it_last = last_states.begin();
    unsigned nth_state(0);
    for( ; it_last != last_states.end(); it_last++ ){

	if( nth_state++ == getSettings()->kalmanMaxNumStatesCutValue() ) break; 

	const kalmanState *state = *it_last;

	if( state->nStubLayers() < getSettings()->minStubLayers() ){
	    break;
	}
	if( last_nstub_layers != 0 && last_nstub_layers != state->nStubLayers() ){
	    if( getSettings()->kalmanSelectMostNumStubState() ) break;
	    best_cands.push_back( cand ); 
	}

	if( cand_chi2 == 0 || cand_chi2 > state->reducedChi2() ){
	    cand = state;
	    cand_chi2 = state->reducedChi2();
	}
	last_nstub_layers = state->nStubLayers();
    }
    if( cand ) best_cands.push_back( cand ); 
    if( getSettings()->kalmanDebugLevel() >= 1 ){
	cout << "------------------------------------" << endl;
	cout << "Best candidates : " << endl; 
	for( auto &s : best_cands ){
	    s->dump( cout, tpa );
	}
	cout << "------------------------------------" << endl;
    }

    //return L1fittedTrk for the selected state.
    if( cand ){

	//candidate dump
	if( getSettings()->kalmanDebugLevel() >= 1 ){

	    cout << "------------------------------------" << endl;
	    if( tpa && tpa->useForAlgEff() ){
		cout << "TP for eff. addr. index : " << tpa << " " << tpa->index() << endl;
	    }
	    cout << "Candidate : " << endl; 
	    cand->dump( cout, tpa, true );
	}

	//fill histograms for the selected state with TP for algEff
	if( getSettings()->kalmanFillInternalHists() ) fillCandHists( *cand, tpa );

	std::map<std::string, double> tp = getTrackParams(cand);
	L1fittedTrack returnTrk(getSettings(), l1track3D, cand->stubs(), tp["qOverPt"], tp["d0"], tp["phi0"], tp["z0"], tp["t"], cand->chi2(), nPar_, iPhiSec, iEtaReg, true);

	if( getSettings()->kalmanDebugLevel() >= 1 ){
	    if( tpa && tpa->useForAlgEff() && returnTrk.getPurity() != 1 ){
		cout << "The candidate is not pure" << endl;
	    }
	    cout << "------------------------------------" << endl;
	}
	return returnTrk;
    }
    else{
	//dump on the missed TP for efficiency calculation.
	if( getSettings()->kalmanDebugLevel() >= 1 ){
	    if( tpa && tpa->useForAlgEff() ){
		cout << "TP for eff. missed addr. index : " << tpa << " " << tpa->index() << endl;
		printStubs( cout, stubs );
		//		printStubAssociatedTPs( cout, stubs );
	    }
	}
	std::vector<const Stub*> nostubs(0);
	L1fittedTrack returnTrk(getSettings(), l1track3D, l1track3D.getStubs(), l1track3D.qOverPt(), 0, l1track3D.phi0(), 0, 0, 9999, nPar_, iPhiSec, iEtaReg, false);
	return returnTrk;
    }
}

unsigned L1KalmanComb::getNextLayer( unsigned state_layer, unsigned next_stub_layer )
{
    if( state_layer < 7 ){

	if(  next_stub_layer < 7 ) return state_layer + 1;
	else if( next_stub_layer < 16 ) return 11; 
	else return 21; 
    }
    else return state_layer + 1; 
}

std::vector<const Stub *> L1KalmanComb::getNextLayerStubs( const kalmanState *state, std::vector<const Stub *> &stubs, unsigned &next_layer )
{
    std::vector<const Stub *> list;

    if( stubs.size() == 0 ) return list;

    unsigned state_layer = state->layerId(); 

    unsigned next_stub_layer = stubs.at(0)->layerId();
    next_layer = getNextLayer( state_layer, next_stub_layer );

    std::vector<const Stub *>::iterator s_it = stubs.begin();
    std::vector<const Stub *>::iterator lastnextstub_it = stubs.end();
    for(; s_it != stubs.end(); s_it++ ){
	const Stub *stub = *s_it;
	if( stub->layerId() == next_layer ){
	    list.push_back(stub);
	    lastnextstub_it = s_it;
	}
	else break;
    }
    if( list.size() != 0 ) stubs.erase( stubs.begin(), lastnextstub_it + 1 );

    return list;
}

std::vector<const kalmanState *> L1KalmanComb::doKF( unsigned nItr, const std::vector<const kalmanState *> &states, std::vector<const Stub *> stubs, const TP *tpa ){
    
    if( getSettings()->kalmanDebugLevel() >= 2 ){
	cout << "----------------------------" << endl;
	cout << "doKF # of iteration = " << nItr << " # of stubs left = " << stubs.size() << " # of the last states = " << states.size() << endl;
	cout << "----------------------------" << endl;
    }

    //finish when there is no more stub or no state
    if( states.size() == 0 ) return states;
    if( stubs.size() == 0 ) return states;

    std::vector<const kalmanState *> active_states;
    std::vector<const kalmanState *> new_states;
    std::vector<unsigned> nvs(3,0);

    std::vector<const kalmanState *>::const_iterator i_state = states.begin();
    for(; i_state != states.end(); i_state++ ){ 

	const kalmanState *the_state = *i_state;
	if( the_state->nVirtualStubs() == getSettings()->kalmanMaxNumVirtualStubs() + 1 && the_state->nStubLayers() >= getSettings()->minStubLayers() ){
	    new_states.push_back( the_state );
	    nvs.at( the_state->nVirtualStubs()-1 )++;
	}
	else{
	    active_states.push_back( the_state );
	}
    }
    if( active_states.size() == 0 ) return new_states;

    //find next layer id from the last state and get next stub list removing those from the original stub list. 
    unsigned next_layer(0);
    std::vector<const Stub *> pre_next_stubs = getNextLayerStubs( active_states.at(0), stubs, next_layer );

    if( getSettings()->kalmanDebugLevel() >= 2 ){
	cout << "# of pre next stubs = " << pre_next_stubs.size() << endl;
    }

    i_state = active_states.begin();
    for(; i_state != active_states.end(); i_state++ ){ 

	if( new_states.size() == getSettings()->kalmanMaxNumStatesCutValue() ) break; 

	const kalmanState *the_state = *i_state;


	//stub cut based on the stub compatibility with the last evaluated state.
	//No cut for the state with less than 3 stubs.
	std::vector<const Stub *> next_stubs;

	if( (int)the_state->nStubLayers() < 3 ){

	    if( pre_next_stubs.size() <= getSettings()->kalmanMaxNumNextStubs() ) 
		next_stubs = pre_next_stubs;
	}
	else{
	    for( unsigned i=0; i < pre_next_stubs.size(); i++ ){

		const Stub * pre_next_stub = pre_next_stubs[i];
//		if( dump_ ) { cout << "stub layerId, sigmaX, sigmaZ = " << pre_next_stub->layerId() << ", " << pre_next_stub->sigmaX() << ", " << pre_next_stub->sigmaZ() << endl; }

		const kalmanState *state = the_state;
		if( nItr == 1 && fitterName_.compare( "KF5ParamsComb" ) == 0 ){
		    const kalmanState *state0 = updateSeedWithStub( *the_state, pre_next_stub );
		    state = state0;
		}

		double e2(0);
		bool pass = validationGate( pre_next_stub, nItr, *state, e2 );
		if( pass ){
		    next_stubs.push_back( pre_next_stub );
		}
		else{
		    if( getSettings()->kalmanDebugLevel() >= 2 ){
			if( tpa && tpa->useForAlgEff() ){
			    set<const TP*> tps = pre_next_stub->assocTPs();
			    if( state->good( tpa ) && tps.find( tpa ) != tps.end() ){
				cout << "A good stub is thrown away." << " e2 = " << e2 << " TPindex = " << tpa->index() << " [eta,phi] = [" << iCurrentPhiSec_ << " , " << iCurrentEtaReg_ << "]" << endl;
				printStub(cout,pre_next_stub);
				validationGate( pre_next_stub, nItr, *state, e2, true );
			    }
			}
		    }
		}
	    }
	}
	if( getSettings()->kalmanDebugLevel() >= 2 ){
	    cout << "# of next stubs = " << next_stubs.size() << endl;
	}

	//stubs are sorted for overlap stub merging.
	if( next_layer < 10 ) 
	    sort( next_stubs.begin(), next_stubs.end(), orderStubsByZ );
	else
	    sort( next_stubs.begin(), next_stubs.end(), orderStubsByR );

	//stub loop
	for( unsigned i=0; i < next_stubs.size() && i < getSettings()->kalmanMaxNumNextStubs() ; i++ ){

	    const Stub * next_stub = next_stubs[i];
//	    if( dump_ ){ cout << "stub (phi,z) = ( " << next_stub->phi() << ", " << next_stub->z() << ")" << endl; } 

	    //For 5 parameter, seed d0 is calculated from stub's bend information.
	    const kalmanState *state = the_state;
	    if( nItr == 1 && fitterName_.compare( "KF5ParamsComb" ) == 0 ){
		const kalmanState *state0 = updateSeedWithStub( *the_state, next_stub );
		state = state0;
	    }

	    //The stubs close to each others are processed one after another as a set of stubs.
	    const kalmanState *new_state = kalmanUpdate( nItr, next_stub, *state, tpa );
	    while( next_stub != next_stubs.back() && isOverlap( next_stub, next_stubs.at(i+1) ) ){
		if( getSettings()->kalmanFillInternalHists() ) 
		    hnmergeStub_->Fill(0);
		next_stub = next_stubs.at(i+1);
		new_state = kalmanUpdate( nItr, next_stub, *new_state, tpa );
		i++;
	    }

	    //state cut
	    if( isGoodState( *new_state ) ){

		nvs.at( new_state->nVirtualStubs() - 1 )++;
		new_states.push_back( new_state );

	    }
	    else{
		if( getSettings()->kalmanDebugLevel() >= 2 ){
		    if( tpa && tpa->useForAlgEff() ){
			if( new_state->good( tpa ) ){
			    cout << "A good state is thrown away." << " rchi2 = " << new_state->reducedChi2() << " TPindex = " << tpa->index() << " [eta,phi] = [" << iCurrentPhiSec_ << " , " << iCurrentEtaReg_ << "]" << endl;
			    new_state->dump(cout, tpa, true );
			}
		    }
		}
	    } 
	}//end of next stub loop

	//A virtual stub is added to all the states with less than the maximum # of virtual stubs. 
	//The state counts the seed as virtual stub. This is not taken into account in the setting kalmanMaxNumVirtualStubs.
	double r = getRofState( next_layer, the_state->xa() );
	if( the_state->nVirtualStubs() - 1 < getSettings()->kalmanMaxNumVirtualStubs() ){
	    const kalmanState *new_state_vs = mkState( nItr, next_layer, r, the_state, the_state->xa(), the_state->pxxa(), 0, the_state->chi2() ); 
	    new_states.push_back( new_state_vs );
	    nvs.at( new_state_vs->nVirtualStubs() - 1 )++;
	}
    }//end of state loop


    //filling the # of states histograms
    if( getSettings()->kalmanFillInternalHists() ) fillEachNumOfVirtualStubStateHists( nItr, nvs.at(0), nvs.at(1), nvs.at(2) );

    if( getSettings()->kalmanDebugLevel() >= 2 ){

	cout << "STATES at Itr = " << nItr << " layerId = " << next_layer << " [" << nvs.at(0) << " " << nvs.at(1) << " " << nvs.at(2) << " " << new_states.size() << "]" << endl;

	if( new_states.size() == 0 && tpa && tpa->useForAlgEff() ){
	    cout << "No state remained for a good track and no more KF at Iteration : " << nItr << endl;
	}
    }

    return doKF( nItr+1, new_states, stubs, tpa );
}

bool L1KalmanComb::validationGate( const Stub *stub, unsigned stub_itr, const kalmanState &state, double &e2, bool debug )const 
{
        
    e2 = 0;
    if( state.stubs().size() < 3 ) return true; 

    std::vector<double> xa     = state.xa();
    TMatrixD            cov_xa = state.pxxa(); 
    if( state.barrel() && !stub->barrel() ){ barrelToEndcap( xa, cov_xa ); }

    std::vector<double> delta = residual(stub, xa );
    TMatrixD f = F( stub, &state );
    TMatrixD h = H(stub);

    TMatrixD pxxm = PxxModel( &state, stub, stub_itr );
    TMatrixD pxxf = HxxH( f, cov_xa ) + pxxm; 
    TMatrixD pddm = PddMeas( stub, &state );
    TMatrixD pddf( pddm.GetNrows(), pddm.GetNcols() );
    pddf = HxxH( h, pxxf ) + pddm; 
    e2 = Chi2( pddf, delta );

    if( debug ){
	cout << "VALIDATION GATE OUTPUT" << endl;
	cout << "State " << endl;
	state.dump(cout);
	cout << "F" << endl;
	f.Print();
	cout << "H" << endl;
	h.Print();
	cout << "xcov prediction" << endl;
	pxxf.Print();
	cout << "mcov" << endl; 
	pddm.Print();
	cout << "HxcovH + pddm " << endl;
	pddf.Print();
	cout << "Residual = ";
	for( unsigned i=0; i < delta.size(); i++ ) cout << delta.at(i) << " ";
	cout << endl;
    }

    return e2 * 0.5 < getSettings()->kalmanValidationGateCutValue();
}

double L1KalmanComb::calcChi2( unsigned itr, const kalmanState &state )const{
    
    if( getSettings()->kalmanDebugLevel() >= 4 ){
	cout << "calcChi2 " << endl;
    }
    double chi2(0), chi2_p(0);
    if( state.last_state() ) chi2 = state.last_state()->chi2();

    const Stub *stub = state.stub();

    if( stub ){

	std::vector<double> delta = residual( stub, state.xa() );
	TMatrixD dcov = PddMeas( stub, &state );
	if( getSettings()->kalmanDebugLevel() >= 4 ){
	    cout << "dcov" << endl;
	    dcov.Print();
	    cout << "xcov" << endl;
	    state.pxxa().Print();
	}
	TMatrixD h = H(stub);
	TMatrixD hxxh = HxxH( h, state.pxxa() );
	if( getSettings()->kalmanDebugLevel() >= 4 ){
	    cout << "h" << endl;
	    h.Print();
	    cout << "hxxh" << endl;
	    hxxh.Print();
	}
	TMatrixD covR = dcov - hxxh;
	if( getSettings()->kalmanDebugLevel() >= 4 ){
	    cout << "covR" << endl;
	    covR.Print();
	    cout << "---" << endl;
	    cout << scientific << "delta = " << delta[0] << ", " << delta[1] << endl;
	}
	chi2_p = Chi2( covR, delta );  
	if( getSettings()->kalmanDebugLevel() >= 3 ){
	    if( chi2_p < 0 ){
		cout << "CALC CHI SQUARE OUTPUT" << endl;
		cout << "WARNING.  CHI SQUARE IS NEGATIVE." << endl;
		cout << "dcov" << endl;
		dcov.Print();
		cout << "xcov" << endl;
		state.pxxa().Print();
		cout << "hxxh" << endl;
		hxxh.Print();
		Chi2( covR, delta, true );
	    }
	}
    }
    chi2 += chi2_p;

    return chi2;
}

const kalmanState *L1KalmanComb::kalmanUpdate( unsigned thisItr, const Stub *stub, const kalmanState &state, const TP *tpa ){

    if( getSettings()->kalmanDebugLevel() >= 3 ){
	cout << "---------------" << endl;
	cout << "kalanUpdate" << endl;
	cout << "---------------" << endl;
	printStub( cout, stub );
    }


    std::vector<double> xa     = state.xa();
    TMatrixD            cov_xa = state.pxxa(); 
    if( state.barrel() && !stub->barrel() ){ 
	barrelToEndcap( xa, cov_xa );
	if( getSettings()->kalmanDebugLevel() >= 3 ){
	    cout << "Previous state changed from Barrel to Endcap parameters" << endl;
	}
    }
    TMatrixD f = F(stub, &state );
    TMatrixD ft(TMatrixD::kTransposed, f );
    if( getSettings()->kalmanDebugLevel() >= 3 ){
	cout << "f" << endl;
	f.Print();
	cout << "ft" << endl;
	ft.Print();
    }

    std::vector<double> fx = Fx( f, xa ); 
    if( getSettings()->kalmanDebugLevel() >= 3 ){
	cout << "fx = ["; 
	for( unsigned i = 0; i < nPar_; i++ ) cout << fx.at(i) << ", ";
	cout << "]" << endl;
    }

    std::vector<double> delta = residual(stub, fx );
    if( getSettings()->kalmanDebugLevel() >= 3 ){
	cout << "delta = " << delta[0] << ", " << delta[1] << endl;
    }

    TMatrixD h = H(stub);
    if( getSettings()->kalmanDebugLevel() >= 3 ){
	cout << "h" << endl;
	h.Print();
    }
    /*
       TMatrixD hth(nPar_, nPar_);
       for( int i=0; i < h.GetNcols(); i++ ){
       for( int j=0; j < h.GetNrows(); j++ ){
       for( int k=0; k < h.GetNcols(); k++ ){
       hth( i, k ) += h( j, i ) * h( j, k );
       }
       }
       }
       if( getSettings()->kalmanDebugLevel() >= 3 ){
       cout << "hth" << endl;
       hth.Print();
       }
       */


    if( getSettings()->kalmanDebugLevel() >= 3 ){
	cout << "previous state covariance" << endl;
	cov_xa.Print();
    }
    TMatrixD pxxm = PxxModel( &state, stub, thisItr );
    if( getSettings()->kalmanDebugLevel() >= 3 ){
	cout << "pxxm" << endl;
	pxxm.Print();
    }

    TMatrixD pxcov = f * cov_xa * ft + pxxm;
    if( getSettings()->kalmanDebugLevel() >= 3 ){
	cout << "pxcov" << endl;
	pxcov.Print();
    }
    TMatrixD dcov = PddMeas( stub, &state );
    if( getSettings()->kalmanDebugLevel() >= 3 ){
	cout << "dcov" << endl;
	dcov.Print();
    }
    TMatrixD k = GetKalmanMatrix( h, pxcov, dcov );  
    if( getSettings()->kalmanDebugLevel() >= 3 ){
	cout << "k" << endl;
	k.Print();
    }
    //   cout << "adjust starts" << endl;
    std::vector<double> new_xa(nPar_);
    TMatrixD new_pxxa;
    GetAdjustedState( k, pxcov, xa, stub, new_xa, new_pxxa );
    if( getSettings()->kalmanDebugLevel() >= 3 ){
	if( nPar_ == 4 )
	    cout << "adjusted x = " << new_xa[0] << ", " << new_xa[1] << ", " << new_xa[2] << ", " << new_xa[3] << endl;
	else if( nPar_ == 5 )
	    cout << "adjusted x = " << new_xa[0] << ", " << new_xa[1] << ", " << new_xa[2] << ", " << new_xa[3] << ", " << new_xa[4] << endl;
	cout << "adjusted covx " << endl;
	new_pxxa.Print();
    }

    std::vector<const Stub*> stubs = state.stubs();
    stubs.push_back(stub);

    const kalmanState *new_state = mkState( thisItr, stub->layerId(), stub->r(), &state, new_xa, new_pxxa, stub, 0 );
    if( getSettings()->kalmanDebugLevel() >= 3 ){
	cout << "new state" << endl;
	new_state->dump( cout, tpa  );
    }

    if( getSettings()->kalmanFillInternalHists() ) fillStepHists( tpa, thisItr, pxcov, pxxm, dcov, dcov + HxxH(h,pxcov), k, new_state );  

    //   cout << "kalanUpdate end" << endl;
    return new_state;
}
void L1KalmanComb::resetStates()
{
    for( unsigned int i=0; i < state_list_.size(); i++ ){

	delete state_list_.at(i);
    }
    state_list_.clear();
}
const kalmanState *L1KalmanComb::mkState( unsigned nIterations, unsigned layerId, double r, const kalmanState *last_state, 
	const std::vector<double> &x, const TMatrixD &pxx, const Stub* stub, double chi2 )
{
    //    cout << "mkState" << endl;
    kalmanState *new_state = new kalmanState( nIterations, layerId, last_state, x, pxx, stub, chi2, this, &getTrackParams );

    if( chi2 == 0 ){
	double new_state_chi2 = calcChi2( nIterations, *new_state ); 
	new_state->setChi2( new_state_chi2 );
    }

    state_list_.push_back( new_state );
    return new_state;
}

std::vector<double> L1KalmanComb::residual(const Stub* stub, const std::vector<double> &x )const{

    std::vector<double> vd = d(stub );
    std::vector<double> hx = Hx( H(stub), x ); 
    std::vector<double> delta(2);
    for( unsigned i=0; i<2; i++ ) delta.at(i) = vd.at(i) - hx.at(i);
    delta.at(0) = wrapRadian(delta.at(0));
    return delta;
}


void L1KalmanComb::bookHists(){

    edm::Service<TFileService> fs_;
    string dirName;
    if( fitterName_.compare("") == 0 ) dirName = "L1KalmanCombInternal";
    else dirName = fitterName_ + "Internal";

    TFileDirectory inputDir = fs_->mkdir(dirName.c_str());

    TString hname;

    hnmergeStub_ = inputDir.make<TH1F>( "hnmergeStub", "# of merged stubs", 1, 0, 1 );
    float nbins(2002);
    for( unsigned i=0; i < nPar_; i++ ){
	hname = Form( "hxt_%d", i );
	hxtMap[hname] = inputDir.make<TH1F>( hname, Form( "; true state values %d", i ), nbins, hxmin[i], hxmax[i] );
	hname = Form( "hx0_%d", i );
	hx0Map[hname] = inputDir.make<TH1F>( hname, Form( "; after HT state values %d", i ), nbins, hxmin[i], hxmax[i] );
	hname = Form( "hxf_%d", i );
	hxfMap[hname] = inputDir.make<TH1F>( hname, Form( "; after KF state values %d", i ), nbins, hxmin[i], hxmax[i] );
    }

    hname = Form( "hrchi2_last_sel_goodtp" );
    hchi2Map[hname] = inputDir.make<TH1F>( hname, Form( "; reduced #Chi^2 after the last iteration for the selected candidate with good TP" ), nbins, hchi2min, hchi2max );
    for( unsigned i=0; i < nPar_; i++ ){
	hname = Form( "hxaxt_last_sel_goodtp_%d", i );
	hxaxtMap[hname] = inputDir.make<TH1F>( hname, Form( "; Estimated x - true x (%d) at last Iteration for the selected candidate with good TP", i ), nbins, hxaxtmin[i], hxaxtmax[i] );
    }
    for( unsigned itr=1; itr <= 10; itr++ ){
	hname = Form( "hchi2_sel_goodtp_itr%d", itr );
	hchi2Map[hname] = inputDir.make<TH1F>( hname, Form( "; #Chi^2 at Iteration %d for the selected candidate with good TP", itr ), nbins, hchi2min, hchi2max );
	for( unsigned i=0; i < nPar_; i++ ){
	    hname = Form( "hxaxt_sel_goodtp_itr%d_%d", itr, i );
	    hxaxtMap[hname] = inputDir.make<TH1F>( hname, Form( "; Estimated x - true x (%d) at Iteration %d for the selected candidate with good TP", i, itr ), nbins, hxaxtmin[i], hxaxtmax[i] );
	    hname = Form( "hkfx_itr%d_%d", itr, i );
	    hkfxMap[hname] = inputDir.make<TH1F>( hname, Form( "; KF x(%d) at Iteration %d for the selected candidate with good TP", i, itr ), nbins, hkfxmin[i], hkfxmax[i] );
	}
    }

    hname = Form( "hrchi2_last" );
    hchi2Map[hname] = inputDir.make<TH1F>( hname, Form( "; reduced #Chi^2 after the last iteration" ), nbins, hchi2min, hchi2max );
    for( unsigned i=0; i < nPar_; i++ ){
	hname = Form( "hxaxt_last_%d", i );
	hxaxtMap[hname] = inputDir.make<TH1F>( hname, Form( "; Estimated x - true x (%d) at last Iteration", i ), nbins, hxaxtmin[i], hxaxtmax[i] );
    }
    for( unsigned itr=1; itr <= 10; itr++ ){
	hname = Form( "hchi2_itr%d", itr );
	hchi2Map[hname] = inputDir.make<TH1F>( hname, Form( "; #Chi^2 at Iteration %d", itr ), nbins, hchi2min, hchi2max );
	for( unsigned i=0; i < nPar_; i++ ){
	    hname = Form( "hxaxt_itr%d_%d", itr, i );
	    hxaxtMap[hname] = inputDir.make<TH1F>( hname, Form( "; Estimated x - true x (%d) at Iteration %d", i, itr ), nbins, hxaxtmin[i], hxaxtmax[i] );
	    for( unsigned j=0; j < nMeas_; j++ ){
		hname = Form( "hk_itr%d_%d_%d", itr, i, j );
		hkMap[hname] = inputDir.make<TH1F>( hname, Form( "; K(%d,%d) at Iteration %d", i, j, itr ), nbins, -1 * hdxmin[i]*hddmin[j], hdxmax[i] * hddmax[j] );
	    }
	}
    }
    for( unsigned i=0; i < nPar_; i++ ){
	for( unsigned j=0; j <= i; j++ ){
	    for( unsigned itr=0; itr<=10; itr++ ){

		hname = Form( "hpxxa_itr%d_%d_%d", itr, i, j );
		hpxxaMap[hname] = inputDir.make<TH1F>( hname, Form( "; state covariance adjusted values on iteration %d (%d,%d)", itr, i, j ), 
			nbins, -1 * hdxmin[i]*hdxmin[j], hdxmax[i]*hdxmax[j] );
		if( itr== 0 ) continue;

		hname = Form( "hpxxf_itr%d_%d_%d", itr, i, j );
		hpxxfMap[hname] = inputDir.make<TH1F>( hname, Form( "; state covariance forcast values on iteration %d (%d,%d)", itr, i, j ), 
			nbins, -1 * hdxmin[i]*hdxmin[j], hdxmax[i]*hdxmax[j] );

		hname = Form( "hpxxModel_itr%d_%d_%d", itr, i, j );
		hPxxModelMap[hname] = inputDir.make<TH1F>( hname, Form( "; state model covariance values on iteration %d (%d,%d)", itr, i, j ), 
			nbins, -1 * hdxModelmin[i]*hdxModelmin[j], hdxModelmax[i]*hdxModelmax[j] );
	    }
	}
    }
    for( unsigned l=0; l<16; l++ ){
	for( unsigned itr=1; itr <= 10; itr++ ){
	    for( unsigned i=0; i < nMeas_; i++ ){
		hname = Form( "hres_itr%d_layer%d_%d", itr, LayerId[l], i );
		hresMap[hname] = inputDir.make<TH1F>( hname, Form( "; residual values on iteration %d (%d) in layer id=%d", itr, LayerId[l], i ), 
			nbins, hresmin[i], hresmax[i] );
		for( unsigned j=0; j <= i; j++ ){
		    hname = Form( "hpddf_itr%d_layer%d_%d_%d", itr, LayerId[l], i, j );
		    hpddfMap[hname] = inputDir.make<TH1F>( hname, Form( "; residual covariance forcast values on iteration %d (%d,%d) in layer id=%d", itr, LayerId[l], i, j ), 
			    nbins, -1 * hddmin[i]*hddmin[j], hddmax[i]*hddmax[j] );
		    hname = Form( "hpddMeas_itr%d_layer%d_%d_%d", itr, LayerId[l], i, j );
		    hPddMeasMap[hname] = inputDir.make<TH1F>( hname, Form( "; measurement error covariance values on iteration %d (%d,%d) in layer id=%d", itr, LayerId[l], i, j ), 
			    nbins, -1 * hddMeasmin[i]*hddMeasmin[i], hddMeasmax[i]*hddMeasmax[j] );
		}
	    }
	}
    }
    hTrackStubEtaLayer = inputDir.make<TH2F>( "hTrackStubEtaLayer", "; stub #eta; stub layer id", 50, -2.5, 2.5, 30, 0, 30 ); 

    for( unsigned itr=1; itr <= 10; itr++ ){
	for( unsigned nv=0; nv <= 2; nv++ ){
	    hname = Form( "hNumStatesItr%d_%dVS", itr, nv );
	    hNumStatesItrMap[hname] = inputDir.make<TH1F>( hname,  Form("; # of states at iteration %d for %d virtual stubs", itr, nv ), 500, -0.5, 499.5 );
	}
    }
}
void L1KalmanComb::fillCandHists( const kalmanState &state, const TP *tpa )
{
    TString hname;

    hname = Form( "hrchi2_last" );
    if( hchi2Map.find(hname) == hchi2Map.end() ){
	cout << hname << " does not exist." << endl;
    }
    else{
	hchi2Map[hname]->Fill( state.reducedChi2() );
    }

    if( !( tpa && tpa->useForAlgEff() ) ) return;

    string type( "sel_goodtp" ); 

    hname = Form( "hrchi2_last_%s", type.c_str() );
    if( hchi2Map.find(hname) == hchi2Map.end() ){
	cout << hname << " does not exist." << endl;
    }
    else{
	hchi2Map[hname]->Fill( state.reducedChi2() );
    }
    const kalmanState *last = &state;
    while( last->nIterations() > 0 ){
	hname = Form( "hchi2_%s_itr%d", type.c_str(), last->nIterations() );
	if( hchi2Map.find( hname ) == hchi2Map.end() ){
	    cout << hname << " does not exist." << endl;
	}
	else hchi2Map[hname]->Fill( last->chi2() );
	last = last->last_state(); 
    }

    std::vector<double> xt(nPar_);
    if( tpa ){

	std::vector<double> xf = state.xa();
	std::map<std::string, double> mxf = getTrackParams( &state );
	std::vector<double> vxf(nPar_);
	vxf[0] = mxf["qOverPt"];
	vxf[1] = mxf["phi0"];
	vxf[2] = mxf["z0"];
	vxf[3] = mxf["t"];
	if( nPar_ == 5 ) vxf[4] = mxf["d0"];

	for( unsigned i=0; i < nPar_; i++ ){
	    hname = Form( "hxf_%d", i );
	    if( hxfMap.find(hname) == hxfMap.end() ){
		cout << hname << " does not exist." << endl;
	    }
	    else hxfMap[hname]->Fill(vxf[i]);
	}
	xt[0] = tpa->qOverPt();
	xt[1] = tpa->phi0();
	xt[2] = tpa->z0();
	xt[3] = tpa->tanLambda();
	if( nPar_ == 5 ) xt[4] = tpa->d0();

	for( unsigned i=0; i < nPar_; i++ ){

	    hname = Form( "hxaxt_last_%s_%d", type.c_str(), i );
	    if( hxaxtMap.find(hname) == hxaxtMap.end() ){
		cout << hname << " does not exist." << endl;
	    }
	    else hxaxtMap[hname]->Fill( vxf[i]-xt[i] );

	    const kalmanState *last = &state;
	    while( last->nIterations() > 0 ){
		std::vector<double> x = last->xa();
		std::map<std::string, double> mx = getTrackParams(last);
		std::vector<double> vx(nPar_);
		vx[0] = mx["qOverPt"];
		vx[1] = mx["phi0"];
		vx[2] = mx["z0"];
		vx[3] = mx["t"];
		if( nPar_ == 5 ) vx[4] = mx["d0"];
		hname = Form( "hxaxt_%s_itr%d_%d", type.c_str(), last->nIterations(), i );
		if( hxaxtMap.find( hname ) == hxaxtMap.end() ){
		    cout << hname << " does not exist." << endl;
		}
		else hxaxtMap[hname]->Fill( vx[i] - xt[i] );
		last = last->last_state(); 
	    }
	}
    }
}
void L1KalmanComb::fillTrackHists( const kalmanState *state, const TP *tpa, std::vector<const Stub *> &stubs )
{
    std::vector<double> x0   = state->xa();
    TMatrixD            pxx0 = state->pxxa();
    //Histogram Fill : seed pxxa 
    for( unsigned i=0; i < nPar_; i++ ){
	for( unsigned j=0; j <= i; j++ ){
	    TString hname = Form( "hpxxa_itr%d_%d_%d", 0, i, j );
	    if( hpxxaMap.find( hname ) == hpxxaMap.end() ){
		cout << hname << " does not exist." << endl;
	    }
	    else hpxxaMap[hname]->Fill( pxx0(i,j) );
	}
    }



    //Histogram Fill : True state from TP 
    if( tpa ){
	std::vector<double> xt(nPar_);
	xt[0] = tpa->qOverPt();
	xt[1] = tpa->phi0();
	xt[2] = tpa->z0();
	xt[3] = tpa->tanLambda();
	if( nPar_ == 5 ) xt[4] = tpa->d0();
	for( unsigned i=0; i < nPar_; i++ ){
	    TString hname = Form( "hxt_%d", i );
	    if( hxtMap.find(hname) == hxtMap.end() ){
		cout << hname << " does not exist." << endl;
	    }
	    else hxtMap[hname]->Fill(xt[i]);
	}
	//Histogram Fill : Seed state 
	std::map<std::string, double> mx0 = getTrackParams( state );
	std::vector<double> vx0(nPar_);
	vx0[0] = mx0["qOverPt"];
	vx0[1] = mx0["phi0"];
	vx0[2] = mx0["z0"];
	vx0[3] = mx0["t"];
	if( nPar_ == 5 ) vx0[4] = mx0["d0"];
	for( unsigned i=0; i < nPar_; i++ ){
	    TString hname = Form( "hx0_%d", i );
	    if( hx0Map.find(hname) == hx0Map.end() ){
		cout << hname << " does not exist." << endl;
	    }
	    else hx0Map[hname]->Fill(vx0[i]);
	}
    }
    //Histogram Fill : Stub layer vs. eta 
    for( unsigned s = 0; s < stubs.size(); s++ ){
	const Stub *stub = stubs[s];
	hTrackStubEtaLayer->Fill( stub->eta(), stub->layerId() );
    }
}
void L1KalmanComb::fillEachNumOfVirtualStubStateHists( unsigned nItr, unsigned nvs0, unsigned nvs1, unsigned nvs2 ) 
{
    TString hname;
    hname = Form( "hNumStatesItr%d_0VS", nItr );
    if( hNumStatesItrMap.find(hname) != hNumStatesItrMap.end() ){
	hNumStatesItrMap[hname]->Fill( nvs0 ); 
    }
    hname = Form( "hNumStatesItr%d_1VS", nItr );
    if( hNumStatesItrMap.find(hname) != hNumStatesItrMap.end() ){
	hNumStatesItrMap[hname]->Fill( nvs1 ); 
    }
    hname = Form( "hNumStatesItr%d_2VS", nItr );
    if( hNumStatesItrMap.find(hname) != hNumStatesItrMap.end() ){
	hNumStatesItrMap[hname]->Fill( nvs2 ); 
    }
}

void L1KalmanComb::fillStepHists( const TP *tpa, unsigned nItr, 
	const TMatrixD &pxxf, const TMatrixD &pxxm, const TMatrixD &pddf,
	const TMatrixD &pddm, const TMatrixD &k, const kalmanState *new_state )
{
    const std::vector<double> &xa = new_state->xa();
    const Stub *stub = new_state->stub();
    const TMatrixD &pxxa = new_state->pxxa();
    double chi2 = new_state->chi2();

    TString hname = Form( "hchi2_itr%d", nItr );
    if( hchi2Map.find(hname) == hchi2Map.end() ){
	cout << hname << " does not exist." << endl;
    }
    else{
	hchi2Map[hname]->Fill( chi2 );
    }

    for( unsigned i=0; i < nPar_; i++ ){
	hname = Form( "hkfx_itr%d_%d", nItr, i );
	if( hkfxMap.find( hname ) == hkfxMap.end() ){
	    cout << hname << " does not exist." << endl;
	}
	hkfxMap[hname]->Fill( xa.at(i) );

	if( tpa ){
	    hname = Form( "hxaxt_itr%d_%d", nItr, i );
	    if( hxaxtMap.find( hname ) == hxaxtMap.end() ){
		cout << hname << " does not exist." << endl;
	    }
	    else{
		std::vector<double> xt(nPar_);
		if( fitterName_.compare( "KF4ParamsCombV2" ) != 0 ){ 
		    xt[0] = tpa->qOverPt();
		    xt[1] = tpa->phi0();
		    xt[2] = tpa->z0();
		    xt[3] = tpa->tanLambda();

		    if( nPar_ == 5 ) xt[4] = tpa->d0();
		}
		else{
		    xt[2] = tpa->qOverPt();
		    xt[3] = xt[2] * tpa->phi0();
		    xt[0] = xt[2] * tpa->tanLambda();
		    xt[1] = tpa->z0() + xt[0] * tpa->phi0();
		}
		std::map<std::string, double> mx = getTrackParams( new_state );
		std::vector<double> vx(nPar_);
		vx[0] = mx["qOverPt"];
		vx[1] = mx["phi0"];
		vx[2] = mx["z0"];
		vx[3] = mx["t"];
		if( nPar_ == 5 ) vx[4] = mx["d0"];

		hxaxtMap[hname]->Fill( vx[i]-xt[i] );
		//if( i==1 && j==1 ) cout << "pxxf(1,1)=" << pxxf(i,j) << endl;
	    }
	}
	for( unsigned j=0; j <= i; j++ ){
	    TString hname = Form( "hpxxf_itr%d_%d_%d", nItr, i, j );
	    if( hpxxfMap.find( hname ) == hpxxfMap.end() ){
		cout << hname << " does not exist." << endl;
	    }
	    else{
		hpxxfMap[hname]->Fill( pxxf(i,j) );
		//if( i==1 && j==1 ) cout << "pxxf(1,1)=" << pxxf(i,j) << endl;
	    }

	    hname = Form( "hpxxa_itr%d_%d_%d", nItr, i, j );
	    if( hpxxaMap.find( hname ) == hpxxaMap.end() ){
		cout << hname << " does not exist." << endl;
	    }
	    else hpxxaMap[hname]->Fill( pxxa(i,j) );

	    hname = Form( "hpxxModel_itr%d_%d_%d", nItr, i, j );
	    if( hPxxModelMap.find( hname ) == hPxxModelMap.end() ){
		cout << hname << " does not exist." << endl;
	    }
	    else hPxxModelMap[hname]->Fill( pxxm(i,j) );
	}
    }
    for( unsigned i=0; i < nPar_; i++ ){
	for( int j=0; j < pddf.GetNrows(); j++ ){
	    TString hname = Form( "hk_itr%d_%d_%d", nItr, i, j );
	    if( hkMap.find( hname ) == hkMap.end() ){
		cout << hname << " does not exist." << endl;
	    }
	    else hkMap[hname]->Fill( k(i,j) );
	}
    }

    for( int i=0; i < pddf.GetNrows(); i++ ){
	for( int j=0; j <= i; j++ ){
	    TString hname = Form( "hpddf_itr%d_layer%d_%d_%d", nItr, stub->layerId(), i, j );
	    if( hpddfMap.find( hname ) == hpddfMap.end() ){
		cout << hname << " does not exist." << endl;
	    }
	    else{
		hpddfMap[hname]->Fill(pddf(i,j));
	    }
	}
    }
    for( int i=0; i < pddm.GetNrows(); i++ ){
	for( int j=0; j < i; j++ ){
	    TString hname = Form( "hpddMeas_itr%d_layer%d_%d_%d", nItr, stub->layerId(), i, j );
	    if( hPddMeasMap.find( hname ) == hPddMeasMap.end() ){
		cout << hname << " does not exist." << endl;
	    }
	    else hPddMeasMap[hname]->Fill( pddm(i,j) );
	}
    }
    std::vector<double> delta_new = residual(stub, xa );
    for( unsigned int i=0; i < delta_new.size(); i++ ){
	TString hname = Form( "hres_itr%d_layer%d_%d", nItr, stub->layerId(), i );
	if( hresMap.find(hname) == hresMap.end() ){
	    cout << hname << " does not exist." << endl;
	}
	else hresMap[hname]->Fill( delta_new[i] );  
    }

}

