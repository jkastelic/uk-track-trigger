///=== This is the Kalman Combinatorial Filter for 4 helix parameters track fit algorithm.
 
///=== Written by: Sioni Summers
 
#include "TMTrackTrigger/TMTrackFinder/interface/KF4ParamsComb.h"
#include "TMTrackTrigger/TMTrackFinder/interface/kalmanState.h"
#define CKF_DEBUG
 
static double wrapRadian( double t ){

    if( t > 0 ){
	while( t > M_PI ) t-= 2*M_PI; 
    }
    else{
	while( t < - M_PI ) t+= 2*M_PI; 
    }
    return t;
}

KF4ParamsComb::KF4ParamsComb(const Settings* settings, const uint nPar, const string &fitterName ) : L1KalmanComb(settings, nPar, fitterName ){

    hdxmin[1] = -1.e-2;
    hdxmax[1] = +1.e-2;
    hdxmin[2] = -1.5;
    hdxmax[2] = +1.5;
    hdxmin[3] = -1.e-0;
    hdxmax[3] = +1.e-0;

    hdxmin[1] = -1.e-2;
    hdxmax[1] = +1.e-2;
    hdxmin[2] = -1.5;
    hdxmax[2] = +1.5;
    hdxmin[3] = -1.e-0;
    hdxmax[3] = +1.e-0;

    hdxModelmin[1] = -1.e-5;
    hdxModelmax[1] = +1.e-5;
    hdxModelmin[1] = -1.e-3;
    hdxModelmax[1] = +1.e-3;
    hdxModelmin[2] = -1.e-2;
    hdxModelmax[2] = +1.e-2;
    hdxModelmin[3] = -1.e-7;
    hdxModelmax[3] = +1.e-7;

    hddmin[1] = -10.;
    hddmax[1] = +10.;

    hddMeasmin[1] = -1.e1;
    hddMeasmax[1] = +1.e1;

    hresmin[1] = -10.;
    hresmax[1] = +10.;


    hxaxtmin[0] = -1.e-3;
    hxaxtmax[0] = +1.e-3;
    hxaxtmin[1] = -1.e-1;
    hxaxtmax[1] = +1.e-1;
    hxaxtmin[2] = -10.;
    hxaxtmax[2] = +10.;
    hxaxtmin[3] = -1.e-0;
    hxaxtmax[3] = +1.e-0;
}

std::map<std::string, double> KF4ParamsComb::getTrackParams(const kalmanState *state )const{

    std::vector<double> x = state->xa();
    std::map<std::string, double> y;
    y["qOverPt"] = x.at(INV2R) / getSettings()->invPtToInvR() * 2.; 
    y["phi0"] = wrapRadian( x.at(PHI0) + sectorPhi() );
    y["z0"] = x.at(Z0);
    y["t"] = x.at(T);
    y["d0"] = 0;
    return y;
}
 
/* The Kalman measurement matrix
 * Here I always measure phi(r), and z(r) */
TMatrixD KF4ParamsComb::H(const Stub* stub)const{
    TMatrixD h(2, 4);
    h(PHI,INV2R) = -stub->r();
    h(PHI,PHI0) = 1;
    h(Z,Z0) = 1;
    h(Z,T) = stub->r();
    return h;
}

TMatrixD KF4ParamsComb::dH(const Stub* stub)const{

    double dr(0);
    if(stub->layerId() > 10){
	dr = stub->sigmaZ();
    }

    TMatrixD h(2, 4);
    h(PHI,INV2R) = -dr;
    h(Z,T) = dr;

    return h;
}
 
/* Seed the state vector */
std::vector<double> KF4ParamsComb::seedx(const L1track3D& l1track3D)const{

    std::vector<double> x(nPar_);
    x[INV2R] = getSettings()->invPtToInvR() * l1track3D.qOverPt()/2;
    x[PHI0]  = wrapRadian( l1track3D.phi0() - sectorPhi() );
    x[Z0]    = l1track3D.z0();
    x[T]     = l1track3D.tanLambda();

    return x;
}
 
/* Seed the covariance matrix
 * Note: 1024 is an arbitrary 'large' value */
TMatrixD KF4ParamsComb::seedP(const L1track3D& l1track3D)const{
    TMatrixD p(4,4);
    /*
    double c = getSettings()->invPtToInvR() / 2; 
    p(0,0) = c * 0.015 * c * 0.015;
    p(1,1) = 0.0050*0.0050;
    p(2,2) = 0.65*0.65;
    p(3,3) = 0.016*0.016; 
    */
    p(INV2R,INV2R) = 1.0e-9;  
    p(PHI0,PHI0) = 1.0e-5; 
    p(Z0,Z0) = 1.0e+1; 
    p(T,T) = 1.0e-2; 
    return p;
}
 
/* The forecast matrix
 * (here equals identity matrix) */
TMatrixD KF4ParamsComb::F(const Stub* stub, const kalmanState *state )const{
    TMatrixD F(4,4);
    for(int n = 0; n < 4; n++)
        F(n, n) = 1;
    return F;
}
 
/* the vector of measurements */
std::vector<double> KF4ParamsComb::d(const Stub* stub )const{
    std::vector<double> meas;
    meas.resize(2);
    meas[0] = wrapRadian( stub->phi() - sectorPhi() );
    meas[1] = stub->z();
    return meas;
}
 
/* Measurement uncertainty */
std::vector<double> KF4ParamsComb::ErrMeas(const Stub* stub, std::vector<double> x )const{

    std::vector<double> meas = d(stub);

    std::vector<double> e(2,0);
    if(stub->layerId() < 10){
	double dphi = stub->sigmaX()/stub->r();
	double dz = stub->sigmaZ();
	e[0] = dphi;
	e[1] = dz;
    }else{
	double phi0 = x.at(PHI0); 

	//calculating the relative errors for l & r ( delta_phi = l / r )
	double delta_phi = meas[0] - phi0;
	double l  = stub->r() * delta_phi;
	double dl = stub->sigmaX();
	double dr = stub->sigmaZ();
	double rdl = dl / l;
	double rdr = dr / stub->r();
	double rdphi = sqrt( rdl * rdl + rdr * rdr );
	double dphi = rdphi * delta_phi;
	double dz = stub->sigmaZ();

	e[0] = dphi;
	e[1] = dz;
    }
    return e;
}
TMatrixD KF4ParamsComb::PddMeas(const Stub* stub, const kalmanState *state )const{

    const std::vector<double> &x = state->xa();
    TMatrixD      xx( 4, 4 ); 
    for(unsigned i=0; i < 4; i++ )
	for(unsigned j=0; j < 4; j++ )
	    xx(i,j) = x[i] * x[j];
    TMatrixD dhcov(2,2);
    if( stub->layerId() > 10 ){
	dhcov = HxxH( dH(stub), xx );
    }

    std::vector<double> e = ErrMeas( stub, state->xa() );
    TMatrixD p(2,2);
    p(PHI, PHI) = e[PHI]*e[PHI];
    p(Z,Z)      = e[Z]*e[Z];

    TMatrixD pddm(2,2);
    pddm = dhcov + p; 
    return pddm;
}

/* State uncertainty */
TMatrixD KF4ParamsComb::PxxModel( const kalmanState *state, const Stub *stub, unsigned stub_itr )const
{

    TMatrixD p(4,4);
    if( getSettings()->kalmanMultiScattFactor() == 0 ) return p;
    p(0,0) = 0.01;
    p(1,1) = 0.01;
    p(2,2) = 0.00001;
    p(3,3) = 0.00001;
    return p;
}

std::string KF4ParamsComb::getParams(){
    return "KF4ParamsComb";
}

/* Determine with a stub belongs/does not belong to a candidate
 * Decision based on hit and state uncertainty */
bool KF4ParamsComb::stubBelongs(const Stub* stub, kalmanState& state, unsigned itr )const{

    return true;

}

bool KF4ParamsComb::isGoodState( const kalmanState &state )const
{
    unsigned nStubs = state.stubs().size();
    bool goodState( true );
    double z0=fabs( state.xa()[Z0] ); 
    if( z0 > 20. ) goodState = false;

    if( nStubs >= 3 && state.reducedChi2() > getSettings()->kalmanStateReducedChi2CutValue() ) goodState=false; 
    if( nStubs >= 3 && state.reducedChi2() < -50. ) goodState=false; 

    return goodState;
}
