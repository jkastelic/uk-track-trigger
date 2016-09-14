#include "TMTrackTrigger/TMTrackFinder/interface/KF4ParamsCombV2.h"
#include "TMTrackTrigger/TMTrackFinder/interface/kalmanState.h"
//#include "TMTrackTrigger/TMTrackFinder/interface/Matrix.h"
#include <TMatrixD.h>
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
KF4ParamsCombV2::KF4ParamsCombV2(const Settings* settings, const std::string &fitterName ) : 
    KF4ParamsComb(settings, 4, fitterName ){

    hkfxmin[0] = -15000.;
    hkfxmax[0] = +15000.;
    hkfxmin[1] = -15000.;
    hkfxmax[1] = +15000.;
    hkfxmin[2] = -15000.;
    hkfxmax[2] = +15000.;
    hkfxmin[3] = -15000.;
    hkfxmax[3] = +15000.;

    hdxmin[0] = -4.e-2;
    hdxmax[0] = +4.e-2;
    hdxmin[1] = -1.001e-1;
    hdxmax[1] = +1.001e-1;
    hdxmin[2] = -5.;
    hdxmax[2] = +5.;
    hdxmin[3] = -4.e-1;
    hdxmax[3] = +4.e-1;

    hddMeasmin[1] = -1.e-1;
    hddMeasmax[1] = +1.e-1;

    hresmin[0] = -100.;
    hresmax[0] = +100.;
    hresmin[1] = -100.;
    hresmax[1] = +100.;

    hxaxtmin[0] = -1000;
    hxaxtmax[0] = +1000;
    hxaxtmin[1] = -1000;
    hxaxtmax[1] = +1000;
    hxaxtmin[2] = -1000;
    hxaxtmax[2] = +1000;
    hxaxtmin[3] = -1000;
    hxaxtmax[3] = +1000;
}

std::string KF4ParamsCombV2::getParams(){
    return "KF4ParamsCombV2";
}

std::map<std::string, double> KF4ParamsCombV2::getTrackParams( const kalmanState *state )const{

    std::vector<double> x = state->xa();

    std::map<std::string, double> z;
    double beta = x.at(BETA);
    double z0p  = x.at(Z0P);
    double R0p  = x.at(R0P);
    double rho0 = x.at(RHO0);

    z["qOverPt"] =  1./( getSettings()->invPtToInvR()  * 0.5 * R0p ); 
    z["phi0"] = wrapRadian( rho0 / R0p + sectorPhi() ); 
    z["z0"] = z0p - beta * z["phi0"]; 
    z["t"] = beta / R0p;
    return z;
}
 
std::vector<double> KF4ParamsCombV2::residual(const Stub* stub, std::vector<double> &x )const
{
    std::vector<double> hx = Hx( H(stub), x ); 
    std::vector<double> vd  = d(stub); 

    std::vector<double> delta(2); 
    for( unsigned i=0; i<2; i++ ){
	delta.at(i) = vd.at(i) - hx.at(i);
    }

    return delta;
}

/* Seed the state vector */
std::vector<double> KF4ParamsCombV2::seedx(const L1track3D& l1track3D)const{
    std::vector<double> x;
    x.resize(4);
    double InvR0 = getSettings()->invPtToInvR() * l1track3D.qOverPt();
    double R0 = 1./InvR0;
    double beta = 2 * R0 * l1track3D.tanLambda();

    x[BETA] = beta; 
    x[Z0P] = l1track3D.z0() + beta * l1track3D.phi0();
    x[R0P] = 2. * R0; 
    x[RHO0] = 2. * R0 * wrapRadian( l1track3D.phi0() - sectorPhi() );
    return x;
}

TMatrixD KF4ParamsCombV2::seedP(const L1track3D& l1track3D)const{
    TMatrixD p(4,4);
    for(int n = 0; n < 4; n++)
    for(int i = 0; i < 4; i++)
	p(n,i) = 100.0;
/*
    p(BETA,BETA) = 1.e-3; 
    p(Z0P,Z0P)   = 1.e-2;
    p(R0P,R0P)   = 5;
    p(RHO0,RHO0) = 1.e-1;
    */

    return p;
}

/* the vector of measurements */
std::vector<double> KF4ParamsCombV2::d(const Stub* stub )const{

    std::vector<double> meas;
    meas.resize(2);
    meas[0] = stub->z();
    meas[1] = stub->r();
    return meas;
}

/* The Kalman measurement matrix
 * Here I always measure phi(r), and z(r) */
TMatrixD KF4ParamsCombV2::H(const Stub* stub)const{
    TMatrixD h(2, 4);
    h(0,0) = -( stub->phi() - sectorPhi() );
    h(0,1) = 1;
    h(1,2) = -( stub->phi() - sectorPhi() );
    h(1,3) = 1;
    return h;
}
TMatrixD KF4ParamsCombV2::dH(const Stub* stub, const kalmanState *state )const{


    double dphi = stub->sigmaX() / stub->r();

    if( !stub->barrel() ){
	std::vector<double> x = state->xa();
	double R0p  = x.at(R0P);
	double rho0 = x.at(RHO0);
	double phi0 = rho0 / R0p; 

	//calculating the relative errors for l & r ( delta_phi = l / r )
	double delta_phi = wrapRadian( ( stub->phi() - sectorPhi() ) - phi0 );
	delta_phi = wrapRadian( delta_phi );
	double l  = stub->r() * delta_phi;
	double dl = stub->sigmaX();
	double dr = stub->sigmaZ();
	double rdl = dl / l;
	double rdr = dr / stub->r();
	double rdphi = sqrt( rdl * rdl + rdr * rdr );
	dphi = rdphi * delta_phi;
    }

    TMatrixD h(2, 4);
    h(0,0) = -dphi;
    h(1,2) = -dphi;
    return h;

}

TMatrixD KF4ParamsCombV2::PxxModel( const kalmanState *state, const Stub *stub, unsigned stub_itr )const{
    //not easy to implement the multiple scattering.
    TMatrixD p(4,4,0);

    return p;
}

/* Measurement uncertainty */
std::vector<double> KF4ParamsCombV2::ErrMeas(const Stub* stub, std::vector<double> x )const{

    std::vector<double> e(2);

    if(stub->layerId() < 10){
	e[0] = stub->sigmaZ();
	e[1] = 1.e-2; 
    }else{
	e[0] = 1.e-2; 
	e[1] = stub->sigmaZ(); 
    }
    return e;
}
TMatrixD KF4ParamsCombV2::PddMeas(const Stub* stub, const kalmanState *state )const
{
	using namespace std;
	
    const std::vector<double> &x = state->xa();
    TMatrixD      xx( 4, 4 ); 
    for(unsigned i=0; i < 4; i++ )
	for(unsigned j=0; j < 4; j++ )
	    xx(i,j) = x[i] * x[j];

    TMatrixD dh = dH(stub, state );
    TMatrixD dhcov(2,2);
    dhcov = HxxH( dh, xx );

    TMatrixD p(2,2);
    std::vector<double> e = ErrMeas( stub, state->xa() );
    p(0,0) = e[0] * e[0];
    p(1,1) = e[1] * e[1];


    if( getSettings()->kalmanDebugLevel() >= 4 ){
	cout << "xx" << endl;
	xx.Print();
	cout << "dh" << endl;
	dh.Print();
	cout << "dhcov" << endl;
	dhcov.Print();
	cout << "mcov" << endl;
	p.Print();
    }

    TMatrixD pddm(2,2);
    pddm = dhcov + p; 
    return pddm;
}


/* Determine with a stub belongs/does not belong to a candidate
 * Decision based on hit and state uncertainty */
bool KF4ParamsCombV2::stubBelongs(const Stub* stub, kalmanState& state, std::vector<double> residual )const{

    std::map<std::string,double> x = getTrackParams( &state ); 
    std::vector<double> e = ErrMeas( stub, state.xa() );

    bool goodMeas( true );
    /*
       if( !( std::abs(residual[0]) < 5 * e[0] && std::abs(residual[1]) < 5 * e[1] ) ){
       goodMeas = false;
       cout << "res: " << residual[0] << "," << residual[1] << ", e : " << e[0] << ", " << e[1] << endl;
       }
       */
    bool beamSpotCompatible = 1;
    /*
       if( fabs( x["z0"] ) > 20){
       beamSpotCompatible = 0;
       }
       */
    return goodMeas && beamSpotCompatible;
}

bool KF4ParamsCombV2::isGoodState( const kalmanState &state )const
{
    unsigned nStubs = state.stubs().size();
    bool goodState( true );
    std::map<std::string, double> x = getTrackParams( &state );
    double z0=fabs( x["z0"] ); 
    if( z0 > 20. ) goodState = false;
    /*
       if( nStubs >= 2 ){
       double z0=fabs( state.xa()[Z0] ); 
       double d0=fabs( state.xa()[D0] ); 
       if( !( z0 < 20. && d0 < 3.5 ) ) goodState = false;
       }
       */

    if( nStubs >= 3 && state.reducedChi2() > getSettings()->kalmanStateReducedChi2CutValue() ) goodState=false; 

    return goodState;

}

