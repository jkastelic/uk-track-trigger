#include "TMTrackTrigger/TMTrackFinder/interface/KF5ParamsComb.h"
#include "TMTrackTrigger/TMTrackFinder/interface/kalmanState.h"
//#include "TMTrackTrigger/TMTrackFinder/interface/Matrix.h"
#include <TMatrixD.h>
#define CKF_DEBUG

static unsigned nlayer_eta[25] = 
{ 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 7, 7, 7,
  7, 7, 7, 7, 6, 6, 6, 6, 6, 6};

static double matx_outer[25] = {
0.16, 0.17, 0.18, 0.19, 0.20, 
0.21, 0.26, 0.22, 0.26, 0.38,
0.41, 0.40, 0.44, 0.50, 0.54,
0.60, 0.44, 0.48, 0.60, 0.68,
0.50, 0.48, 0.64, 0.39, 0.20
};

static double matx_inner[25] = {
0.14, 0.1, 0.1, 0.1, 0.1, 
0.1, 0.1, 0.1, 0.1, 0.1, 
0.12, 0.0, 0.1, 0.1, 0.15,
0.20, 0.25, 0.25, 0.3, 0.3,
0.35, 0.40, 0.40, 0.6, 0.6
};

static double wrapRadian( double t ){

    if( t > 0 ){
	while( t > M_PI ) t-= 2*M_PI; 
    }
    else{
	while( t < - M_PI ) t+= 2*M_PI; 
    }
    return t;
}
KF5ParamsComb::KF5ParamsComb(const Settings* settings, const string &fitterName ) : 
    L1KalmanComb(settings, 5, fitterName ){

    hxmin[0] = -0.05;
    hxmax[0] = +0.05;
    hxmin[1] = -3.2;
    hxmax[1] = +3.2;
    hxmin[2] = -20;
    hxmax[2] = +20;
    hxmin[3] = -10;
    hxmax[3] = +10;
    hxmin[4] = -3.5;
    hxmax[4] = +3.5;

    hdxmin[0] = -4.004e-5;
    hdxmax[0] = +4.004e-5;
    hdxmin[1] = -4.004e-3;
    hdxmax[1] = +4.004e-3;
    hdxmin[2] = -4.001e+0;
    hdxmax[2] = +4.001e+0;
    hdxmin[3] = -1.001e-1;
    hdxmax[3] = +1.001e-1;
    hdxmin[4] = -1.001;
    hdxmax[4] = +1.001;

    hdxModelmin[1] = -1.e-5;
    hdxModelmax[1] = +1.e-5;
    hdxModelmin[1] = -1.e-3;
    hdxModelmax[1] = +1.e-3;
    hdxModelmin[2] = -1.e-2;
    hdxModelmax[2] = +1.e-2;
    hdxModelmin[3] = -1.e-7;
    hdxModelmax[3] = +1.e-7;
    hdxModelmin[4] = -1.e-7;
    hdxModelmax[4] = +1.e-7;

    hddmin[1] = -10.;
    hddmax[1] = +10.;

    hddMeasmin[1] = -1.e1;
    hddMeasmax[1] = +1.e1;

    hresmin[1] = -10.;
    hresmax[1] = +10.;


    hxaxtmin[0] = -1.001e-4;
    hxaxtmax[0] = +1.001e-4;
    hxaxtmin[1] = -1.001e-2;
    hxaxtmax[1] = +1.001e-2;
    hxaxtmin[2] = -1.001e+1;
    hxaxtmax[2] = +1.001e+1;
    hxaxtmin[3] = -1.001e-2;
    hxaxtmax[3] = +1.001e-2;
    hxaxtmin[4] = -1.001;
    hxaxtmax[4] = +1.001;
}

std::string KF5ParamsComb::getParams(){
    return "KF5ParamsComb";
}

std::map<std::string, double> KF5ParamsComb::getTrackParams(const kalmanState *state )const{

    std::vector<double> x = state->xa();
    std::map<std::string, double> y;
    y["qOverPt"] = x.at(INV2R) / getSettings()->invPtToInvR() * 2.; 
    y["phi0"] = wrapRadian( x.at(PHI0) + sectorPhi() );
    y["z0"] = x.at(Z0);
    y["t"] = x.at(T);
    y["d0"] = x.at(D0);
    return y;
}

std::vector<double> KF5ParamsComb::residual(const Stub* stub, const std::vector<double> &x )const
{
    std::vector<double> hx = Hx( H(stub), x ); 
    std::vector<double> vd  = d(stub ); 
    //    cout << "residual " << endl;
    //    cout << "vd = " << vd[0] << ", " << vd[1] << endl;
    //H(stub).Print();
    //    cout << "hx = " << hx[0] << ", " << hx[1] << endl;

    std::vector<double> delta(2); 
    for( unsigned i=0; i<2; i++ ){
	delta.at(i) = vd.at(i) - hx.at(i);
    }
    delta.at(0) = wrapRadian(delta.at(0));
    //   cout << "delta = " << delta[0] << ", " << delta[1] << endl;

    return delta;
}

/* Seed the state vector */
std::vector<double> KF5ParamsComb::seedx(const L1track3D& l1track3D)const{
    std::vector<double> x;
    x[INV2R] = getSettings()->invPtToInvR() * l1track3D.qOverPt()/2;
    x[PHI0]  = wrapRadian( l1track3D.phi0() - sectorPhi() );
    x[Z0]    = l1track3D.z0();
    x[T]     = l1track3D.tanLambda();
    x[D0]    = 0;
    return x;
}

/* Seed the covariance matrix
 * Note: 1024 is an arbitrary 'large' value */
TMatrixD KF5ParamsComb::seedP(const L1track3D& l1track3D)const{
    TMatrixD p(5,5);
    for(int n = 0; n < 5; n++)
	p(n,n) = 0.0;
    //    double c0(1.e6);
    //   double c1(1.e6);
    /*
       double c0(1.0);
       double c1(1.0);
       double c = getSettings()->invPtToInvR() / 2; 
       p(INV2R,INV2R) = c0*c * 0.015 * c * 0.015;
       p(PHI0,PHI0) = c0*0.005*0.005;
       p(Z0,Z0) = c1*0.65*0.65;
       p(T,T) = c1*0.016*0.016; 
    //    p(D0,D0) = c1*3.5*3.5; 
    p(D0,D0) = 1.*1.; 
    double hz(15.);
    if( !getSettings()->useSeedFilter() ) p(Z0,Z0) = hz * hz; 
    */

    p(INV2R,INV2R) = 1.0e-9;  
    p(PHI0,PHI0) = 1.0e-5; 
    p(Z0,Z0) = 1.0e+1; 
    p(T,T) = 1.0e-2; 
    p(D0,D0) = 1.0; 
    return p;
}
/* The forecast matrix
 * (here equals identity matrix) */
TMatrixD KF5ParamsComb::F(const Stub* stub, const kalmanState *state )const{
    TMatrixD F(nPar_,nPar_);
    for(unsigned n = 0; n < nPar_; n++)
	F(n, n) = 1;
    return F;
}

/* the vector of measurements */
std::vector<double> KF5ParamsComb::d(const Stub* stub )const{

    std::vector<double> meas;
    meas.resize(2);
    meas[PHI] = wrapRadian( stub->phi() - sectorPhi() );
    meas[Z] = stub->z();
    return meas;
}

/* The Kalman measurement matrix
 * Here I always measure phi(r), and z(r) */
TMatrixD KF5ParamsComb::H(const Stub* stub)const{
    TMatrixD h(2, 5);
    h(PHI,INV2R) = -stub->r();
    h(PHI,PHI0) = 1;
    if( stub->r() == 0 ) h(PHI,D0) = 99999.;
    else h(PHI,D0) = -1./stub->r();
    // else h(0,4) = 1./stub->r();
    h(Z,Z0) = 1;
    h(Z,T) = stub->r();

    return h;
}
TMatrixD KF5ParamsComb::dH(const Stub* stub)const{

    double dr(0);
    if(stub->layerId() > 10){
	dr = stub->sigmaZ();
	//	dr = stub->rErr();
    }

    TMatrixD h(2, 5);
    h(PHI,INV2R) = -dr;
    if( stub->r() == 0 ) h(PHI,D0) = 99999.;
    else h(PHI,D0) = 1./(stub->r()*stub->r()) * dr;
    h(Z,T) = dr;

    return h;
}

/* State uncertainty */
TMatrixD KF5ParamsComb::PxxModel( const kalmanState *state, const Stub *stub, unsigned stub_itr )const
{

    const std::vector<double> &xa = state->xa();
    unsigned last_update_itr(0);
    double last_update_r(0);
    const kalmanState *last_update_state = state->last_update_state(); 
    if( last_update_state ) {
	last_update_itr = last_update_state->nIterations();
	last_update_r = last_update_state->r();
    }
    double eta = stub->eta();
    unsigned n_state_updates = state->nStubLayers();

    TMatrixD p(nPar_,nPar_);




    if( getSettings()->kalmanMultiScattFactor() == 0 ) return p;

    //multiple scattering

    TMatrixD plambda(nPar_,nPar_);
    TMatrixD pphi(nPar_,nPar_);

    double r = last_update_r;
    double dtheta0;

    double dl_inner(0), dl_outer(0);
    unsigned i_eta = abs( eta / 0.1 );
    if( i_eta > 24 ) i_eta = 24;


    dl_inner = matx_inner[i_eta];

    double dl1 = matx_outer[i_eta] / ( nlayer_eta[i_eta] + 1 );
    double w(0.5);

    /*
       if( stub_itr <= 4 )
       dl_outer = ( stub_itr - 1 ) * ( 1 + w ) * dl1; 
       else
       dl_outer = 3 * ( 1. + w ) * dl1 + ( stub_itr - 4 ) * ( 1. - w ) * dl1;
       */



    if( last_update_itr == 0  && n_state_updates != 0 ){
	dl_inner = matx_inner[i_eta];
	r = 2.5;
	if( stub_itr <= 4 )
	    dl_outer = ( stub_itr - 1 ) * ( 1. + w ) * dl1; 
	else
	    dl_outer = 3 * ( 1. + w ) * dl1 + ( stub_itr - 4 ) * ( 1. - w ) * dl1;
    }

    else if( n_state_updates == 3 ){
	dl_inner = matx_inner[i_eta];
	r = 2.5;
	if( stub_itr <= 4 )
	    dl_outer = ( stub_itr - 1 ) * ( 1. + w ) * dl1; 
	else
	    dl_outer = 3 * ( 1. + w ) * dl1 + ( stub_itr - 4 ) * ( 1. - w ) * dl1;
    }
    else{
	double cp(1), cs(1);
	if( stub_itr <= 4 ){
	    double c = sqrt( ( stub_itr - 1 ) * ( stub_itr - 1 ) - ( last_update_itr - 1 ) * ( last_update_itr - 1 ) );
	    //double c = stub_itr - last_update_itr;
	    dl_outer = c * ( 1. + w ) * dl1;
	}
	else{
	    if( last_update_itr <= 3 ){
		double cp = sqrt( 9 - ( last_update_itr - 1 ) * ( last_update_itr - 1 ) );
		double cs = sqrt( ( stub_itr - 1 ) * ( stub_itr - 1 ) - 9 );
		//double cp = 4 - last_update_itr;
		//double cs = stub_itr - 4;

		dl_outer = cp * ( 1. + w ) * dl1 + cs * ( 1. - w ) * dl1;
	    }
	    else{
		double c = sqrt( ( stub_itr - 1 ) * ( stub_itr - 1 ) - ( last_update_itr - 1 ) * ( last_update_itr - 1 ) );
		//double c = stub_itr - last_update_itr;
		dl_outer = c * ( 1. - w ) * dl1;
	    }
	}
    }

    double dl = dl_inner + dl_outer;

    std::map<std::string, double> y = getTrackParams(state);
    dtheta0 = 1./sqrt(3) * 0.0136 * (2.*fabs(y["2rInv"]) ) / getSettings()->invPtToInvR() * sqrt(dl)*( 1+0.038*log(dl) ); 
    dtheta0 *= getSettings()->kalmanMultiScattFactor();

    //lambda
    double dlambda = - dtheta0;
    std::vector<double> e_lambda(nPar_, 0);
    e_lambda[INV2R] = y["2rInv"] * y["t"] * dlambda; 
    e_lambda[Z0] = -1 * r * ( 1 + y["t"] * y["t"] ) * dlambda;
    e_lambda[T] = ( 1 + y["t"] * y["t"] ) * dlambda;

    for( unsigned i = 0; i < nPar_; i++ ){
	for( unsigned j = 0; j < nPar_; j++ ){
	    plambda(i,j) = e_lambda[i] * e_lambda[j];  
	}
    }
    //phi
    std::vector<double> e_phi(nPar_, 0);
    e_phi[PHI0] = dtheta0;
    e_phi[D0] = -1. * r * dtheta0;
    //    e_phi[4] = r * dtheta0;
    for( unsigned i = 0; i < nPar_; i++ ){
	for( unsigned j = 0; j < nPar_; j++ ){
	    pphi(i,j) = e_phi[i] * e_phi[j];  
	}
    }

    p = plambda + pphi;

    return p;
}

/* Measurement uncertainty */
std::vector<double> KF5ParamsComb::ErrMeas(const Stub* stub, std::vector<double> x )const{

    std::vector<double> meas = d(stub);

    std::vector<double> e(2,0);
    if(stub->layerId() < 10){
	double dphi = stub->sigmaX()/stub->r();
	double dz = stub->sigmaZ();
	e[PHI] = dphi;
	e[Z] = dz;
    }else{
	double delta_phi = meas[0] - x[PHI0];
	double dr = stub->sigmaZ() *sqrt(12) * 0.5;
	double dz   = dr * std::fabs(x[T]);
	double rdl = stub->sigmaX() / ( stub->r() * delta_phi );
	double rdr = dr / stub->r();
	double rdphi = sqrt( rdl * rdl + rdr * rdr );
	double dphi = rdphi * delta_phi;
	e[PHI] = dphi;
	e[Z] = dz;
    }
    return e;
}

TMatrixD KF5ParamsComb::PddMeas(const Stub* stub, const kalmanState *state )const{

    const std::vector<double> &x = state->xa();
    TMatrixD      xx( 5, 5 ); 
    for(unsigned i=0; i < 5; i++ )
	for(unsigned j=0; j < 5; j++ )
	    xx(i,j) = x[i] * x[j];
    TMatrixD dhcov(2,2);
    if( stub->layerId() > 10 ){
	dhcov = HxxH( dH(stub), xx );
    }
    //dhcov.Print();

    std::map<std::string, double> y = getTrackParams(state);
    TMatrixD p(2,2);
    if(stub->layerId() < 10){
	double dphi = stub->sigmaX()/stub->r();
	double dz = stub->sigmaZ();
	p(PHI,PHI) = dphi * dphi;
	p(Z,Z) = dz * dz;
    }else{
	double dphi = stub->sigmaX()/stub->r();
	p(PHI,PHI) = dphi * dphi;
    }
    TMatrixD pddm(2,2);
    pddm = dhcov + p; 
    return pddm;

    /*
       const std::vector<double> &x = state->xa();
       const TMatrixD      &xcov = state->pxxa();

       std::map<std::string, double> y = vecToMap(x);
       TMatrixD p(2,2,0);
       if(stub->layerId() < 10){
       double dphi = stub->sigmaX()/stub->r();
       double dz = stub->sigmaZ();
       p(PHI,PHI) = dphi * dphi;
       p(Z,Z) = dz * dz;
       }else{
       double dr = stub->sigmaZ();
       double dphi_dr = - dr * y["2rInv"] + y["d0"]/(stub->r()*stub->r()) * dr ;
       double dz_dr   = dr * y["t"];
       double dphidz_dr = dphi_dr * dz_dr;
       double dphi = stub->sigmaX()/stub->r();
       p(PHI,PHI) = dphi_dr * dphi_dr + dphi * dphi + stub->r() * stub->r() * xcov(0,0) + xcov(1,1) + 1./(stub->r()*stub->r()) * xcov(4,4);  
       p(Z,Z) = dz_dr * dz_dr + stub->r()*stub->r() *xcov(3,3); 
       p(PHI,Z) = dphidz_dr; 
       p(Z,PHI) = p(PHI,Z);
       }
       */
    return p;
}


/* Determine with a stub belongs/does not belong to a candidate
 * Decision based on hit and state uncertainty */
bool KF5ParamsComb::stubBelongs(const Stub* stub, kalmanState& state, unsigned nItr )const{

    return true;
}

const kalmanState *KF5ParamsComb::updateSeedWithStub( const kalmanState &state, const Stub *stub )
{
    std::vector<double> xa = state.xa();
    TMatrixD      pxxa = state.pxxa();
    //    xa[4] = -1 * stub->dphi() * stub->r() * xa[0]/xa[0]; 

    double c0(100);
    double dd0 = 2 * stub->dphiRes() * stub->r();
    pxxa(4,4) = c0 * dd0 * dd0; 
    const kalmanState *new_state = mkState( state.nIterations(), state.layerId(), state.r(), state.last_state(), xa, pxxa, state.stub(), state.chi2() );
    return new_state;
}

bool KF5ParamsComb::isGoodState( const kalmanState &state )const
{
    unsigned nStubs = state.stubs().size();
    bool goodState( true );
    double z0=fabs( state.xa()[Z0] ); 
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

double KF5ParamsComb::getRofState( unsigned layerId, const vector<double> &xa )const
{
    double r(0), z(0);

    double tlambda = xa[T];
    double z0 = xa[Z0]; 

    switch( layerId ){

	case 1:
	    return 23.0000;
	case 2:
	    return 35.7168;
	case 3:
	    return 50.8000;
	case 4:
	    return 68.6000;
	case 5:
	    return 88.7901;
	case 6:
	    return 108.0000;
	case 11:
	case 21:
	    z = 134.9445;
	    break;
	case 12:
	case 22:
	    z = 159.7452;
	    break;
	case 13:
	case 23:
	    z = 189.1039;
	    break;
	case 14:
	case 24:
	    z = 223.8583;
	    break;
	case 15:
	case 25:
	    z = 265.0000;
	    break;
    }
    if( 10 < layerId && layerId < 20 ) z *= -1; 

    r = ( z - z0 ) / tlambda; 
    return r;
}

