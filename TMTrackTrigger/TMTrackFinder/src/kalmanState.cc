#include "TMTrackTrigger/TMTrackFinder/interface/kalmanState.h"
//#include "TMTrackTrigger/TMTrackFinder/interface/Matrix.h"
#include <TMatrixD.h>

kalmanState::kalmanState(): nIterations_(0), layerId_(0), xa_(0), pxxa_(), stub_(0), chi2_(0), n_stubs_(0), n_virtual_stubs_(1), fitter_(0), fXtoTrackParams_(0), barrel_(true){
}

kalmanState::kalmanState( unsigned nIterations, unsigned layerId, const kalmanState *last_state, const std::vector<double> &x, const TMatrixD &pxx, const Stub* stub, double chi2,
	L1KalmanComb *fitter, GET_TRACK_PARAMS f ){

    nIterations_ = nIterations;
    layerId_ = layerId;
    last_state_ = last_state;
    xa_ = x;
    pxxa_.Clear();
    pxxa_.ResizeTo( pxx.GetNrows(), pxx.GetNcols() );
    pxxa_ = pxx;
    stub_ = stub;
    chi2_ = chi2;

    const kalmanState *state = this;
    n_stubs_ = 0;
    n_virtual_stubs_ = 0;
    barrel_ = true;
    r_ = 0;
    z_ = 0;
    bool last_stub(false);

    while( state ){
	if( state->stub() ){
	    n_stubs_ ++; 
	    if( !state->stub()->barrel() ) barrel_ = false;
	    if( !last_stub ){
		r_ = state->stub()->r();
		z_ = state->stub()->z();
		last_stub = true;
	    }
	}
	else n_virtual_stubs_++;
	state = state->last_state();
    }
    n_stub_layers_ = nIterations_ + 1 - n_virtual_stubs_;
    fitter_ = fitter;
    fXtoTrackParams_ = f;


}

kalmanState::kalmanState(const kalmanState &p){

    nIterations_ = p.nIterations();
    layerId_ = p.layerId();
    r_ = p.r();
    z_ = p.z();
    last_state_ = p.last_state();
    xa_ = p.xa();
    pxxa_ = p.pxxa();
    stub_ = p.stub();
    chi2_ = p.chi2();
    n_stubs_         = p.nStubs();
    n_virtual_stubs_ = p.nVirtualStubs();
    n_stub_layers_   = p.nStubLayers();
    fitter_ = p.fitter();
    fXtoTrackParams_ = p.fXtoTrackParams();
    barrel_ = p.barrel();
}

kalmanState & kalmanState::operator=( const kalmanState &other )
{
    if (&other == this)
	return *this;

    nIterations_ = other.nIterations();
    layerId_ = other.layerId();
    r_ = other.r();
    z_ = other.z();
    last_state_ = other.last_state();
    xa_ = other.xa();
    pxxa_ = other.pxxa();
    stub_ = other.stub();
    chi2_ = other.chi2();
    n_stubs_ = other.nStubs();
    n_virtual_stubs_ = other.nVirtualStubs();
    n_stub_layers_ = other.nStubLayers();
    fitter_ = other.fitter();
    fXtoTrackParams_ = other.fXtoTrackParams();
    barrel_ = other.barrel();
    return *this;
}

bool kalmanState::good( const TP *tp )const{

    const kalmanState *state = this;
    while( state ){
	const Stub *stub = state->stub();
	if( stub ){
	    set<const TP*> tps = stub->assocTPs();

	    if( tps.find(tp) == tps.end() ) return false; 
	}
	state = state->last_state();
    }
    return true;
}

double kalmanState::reducedChi2() const
{ 
    if( 2 * n_stubs_ - xa_.size() > 0 ) return chi2_ / ( 2 * n_stubs_ - xa_.size() ); 
    else return 0; 
} 

const kalmanState *kalmanState::last_update_state()const
{
    const kalmanState *state = this;
    while( state ){
	if( state->stub() ) return state;
	state = state->last_state();
    }
    return 0;
}
std::vector<const Stub *> kalmanState::stubs()const
{
    std::vector<const Stub *> stubs;

    const kalmanState *state = this;
    while( state ){
	if( state->stub() ) stubs.push_back( state->stub() ); 
	state = state->last_state();
    }
    return stubs;
}

bool kalmanState::order(const kalmanState *left, const kalmanState *right){ return (left->nStubs() > right->nStubs()); }
bool kalmanState::orderReducedChi2(const kalmanState *left, const kalmanState *right){ 
    return ( left->reducedChi2() < right->reducedChi2() );
}

void kalmanState::dump( ostream &os, const TP *tp, bool all )const
{
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
    std::map<std::string, double> y = fXtoTrackParams_( fitter_, this );

    os << "kalmanState : ";
    os << "# of iterations = " << nIterations_ << ", ";
    os << "layerId = " << layerId_ << ", ";
    os << "barrel = " << barrel_ << ", ";
    os << "r = " << r_ << ", "; 
    os << "z = " << z_ << ", ";
    for( auto pair : y ){
	os << pair.first << ":" << y[pair.first] << " "; 
    }
    os << endl;
    os << "xa = ( ";
    for( unsigned i=0; i<xa_.size()-1; i++ ) os << xa_[i] << ", ";
    os << xa_.back() << " )" << endl;

    os << "xcov" << endl;
    pxxa_.Print(); 
    os << " chi2 = " << chi2_ << ", "; 
    os << " # of stubs = " << nStubs() << ", "; 
    os << " # of stublayers = " << nStubLayers() << endl;
    for( auto &stub : stubs() ){
	os << "              stub ";
	os << "[" << stub << "] "; 
	os << "index : " << stub->index() << " ";
	os << "layerId : " << stub->layerId() << " ";
	os  << "[r,phi,z] = ";
	os << "[" << stub->r() << ", " << stub->phi() << ", " << stub->z() << "] ";
	os << " assoc TP indices = [ "; 
	std::set<const TP*> tps = stub->assocTPs();
	for( auto tp : tps ) os << tp->index() << " "; 
	os << "] ";
	os << endl;
    }
    if( tp ){
	os << "\tTP index = " << tp->index() << " useForAlgEff = " << useForAlgEff << " ";
	os << "rel. residual ";
	for( auto pair : tp_x ){
	    os << pair.first << ":" << ( y[pair.first] - pair.second ) / pair.second << " "; 
	}
    }
    else{
	os << "\tTP index = "; 
    }
    os << endl;

    if( stub_ ){
	os << "\tstub [r,phi,z] = ";
	os << "[" << stub_->r() << ", " << stub_->phi() << ", " << stub_->z() << "] ";
	os << " assoc TP indices = [ "; 
	std::set<const TP*> tps = stub_->assocTPs();
	for( auto tp : tps ) os << tp->index() << " "; 
	os << "] ";
    }
    else{
	os << "\tvirtual stub";
    }
    os << endl;

    if( all ){
	const kalmanState *state = last_state();
	if( state ){
	    state->dump( os, tp, all );
	    state = state->last_state();
	}
	else return;
    }
}



