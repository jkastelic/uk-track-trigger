///=== This is the linearized chi2 track fit algorithm.
///=== It is based on the tracklet group's track fit software
///=== https://github.com/EmanuelPerez/cmssw/blob/TTI_62X_TrackTriggerObjects/SimDataFormats/SLHC/interface/L1TTrack.hh

///=== Written by: Alexander D. Morton

#include "TMTrackTrigger/TMTrackFinder/interface/TrackFitLinearAlgo.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Utility.h"
 
//=== Set configuration parameters.
 
TrackFitLinearAlgo::TrackFitLinearAlgo( const Settings* settings, const uint nPar ) : TrackFitGeneric(settings), settings_(settings), print(settings_->debug() ) {
  nPar_ = nPar; 
}


//=== Initialize constants at start of run.

void TrackFitLinearAlgo::initRun() {
  invPtToInvR_ = settings_->invPtToInvR();  // Conversion factor from 1/Pt to 1/radius_of_curvature.
  //--- These parameters configure how the fit is done.
 
  numFittingIterations_ = settings_->numTrackFitIterations(); // Required number of iterations of fit algo.
  killTrackFitWorstHit_ = settings_->killTrackFitWorstHit(); // Optionally kill hit with worst residual.
  generalResidualCut_   = settings_->generalResidualCut(); // The cut used to remove bad stubs (if nStubs > minLayers)
  killingResidualCut_   = settings_->killingResidualCut(); // The cut used to kill off tracks entirely

  //--- These two parameters are used to check if after the fit, there are still enough stubs on the track
 
  // for it to be valid (in case the fit removed stubs with bad residuals).
  minStubLayers_ = settings_->minStubLayers(); // Min. number of layers that must have stubs for track to be declared found.
  minPtToReduceLayers_ = settings_->minPtToReduceLayers(); // Reduce MinStubLayers by 1 for tracks with pt above this cut.
  std::stringstream lSStr;
  lSStr << "linearised" << nPar_ << "Fit_" << numFittingIterations_ << "Iterations_KillWorseHits" << killTrackFitWorstHit_;
  configParameters_ = (lSStr.str());
}

std::string TrackFitLinearAlgo::getParams() {
	return configParameters_;
}
 
//=== Fit a track candidate obtained from the Hough Transform.
//=== Specify which phi sector and eta region it is in.
L1fittedTrack TrackFitLinearAlgo::fit(const L1track3D& l1track3D, unsigned int iPhiSec, unsigned int iEtaReg) {
 
  rinv_= invPtToInvR_*l1track3D.qOverPt();
  phi0_= l1track3D.phi0();
  z0_ = l1track3D.z0();
  t_ = l1track3D.tanLambda();
  d0_ = 0.0;
  stubs_ = l1track3D.getStubs(); 
 
  // Cheat by using MC truth to initialize helix parameters. Useful to check if convergence is the problem.
  
/*  const TP* tp = l1track3D.getMatchedTP();
  if (tp != nullptr) {
    rinv_= invPtToInvR_*tp->qOverPt();
    phi0_= tp->phi0();
    z0_ = tp->z0();
    t_ = tp->tanLambda();
  }
*/
 
  if (stubs_.size() > __MAX_STUBS_PER_TRK__) {
    std::cout << "TrackFitLinearAlgo WARNING: too many stubs on track - skipping it! " << stubs_.size() << std::endl;
    // Return fitted track created directly from input track candidate.
    std::cout << "Too many stubs." << std::endl;
    return L1fittedTrack(settings_, l1track3D, stubs_, l1track3D.qOverPt(), 0., l1track3D.phi0(), l1track3D.z0(), l1track3D.tanLambda(), 999999., 4, iPhiSec, iEtaReg);
  }
 
  float largestresid;
  int ilargestresid;
  largestresid=-1.0;
  ilargestresid=-1;

  calculateDerivatives( false ); // Default = false
  linearTrackFit( false ); // Default = false
  residuals( largestresid,ilargestresid );

  // If the option to kill the hit with the worse residual is enabled, this is down if the following conditions are met:

  for (int i=1;i<numFittingIterations_+1;++i) {
    if (i>1) {
      if (print) std::cout << __LINE__ << " - killTrackFitWorstHit_= " << killTrackFitWorstHit_ << "/largestresid = " << largestresid << std::endl;
      if (killTrackFitWorstHit_) {
        if ( largestresid > min(killingResidualCut_, generalResidualCut_) ) {
          const vector<const Stub*> stubs_backup(stubs_);
          stubs_.erase(stubs_.begin()+ilargestresid);    
          if (largestresid > killingResidualCut_ || Utility::countLayers( settings_, stubs_ ) >= minStubLayers_) {
            if (print) std::cout << "Killed stub " << ilargestresid << "." << std::endl;
          } else {
            // Don't delete worst stub, as it would kill the track.
            stubs_ = stubs_backup;
          }
        }       
      }
      rinv_=rinvfit4par_;
      phi0_=phi0fit4par_;
      z0_=z0fit4par_;
      d0_=d0fit_;
      t_=tfit4par_;
      
      //      if (print && (i==numFittingIterations_+1)) std::cout << "Last fit." << std::endl;
      
      calculateDerivatives( false ); // Default = false
      linearTrackFit( false ); // Default = false
      residuals(largestresid,ilargestresid);
    }
  }  
    
  // Cheat by using MC truth to initialize helix parameters. Useful to check if convergence is the problem.
  /*
  if (tp != nullptr) {
    rinv_= invPtToInvR_*tp->qOverPt();
    phi0_= tp->phi0();
    z0_ = tp->z0();
    t_ = tp->tanLambda();
  }
  */

  //  if (print) std::cout << "5 fit." << std::endl;
  if (nPar_ == 5){
    calculateDerivatives( true );
    linearTrackFit( true );
  }
 
  // Note if 4 and 5 parameter helix fits give valid fitted tracks.
  // They may not, if the fits removed too many stubs (with bad residuals) from the track.
  // Currently, since both 4 and 5 parameter fits yield the same stubs, they are both valid, or both invalid.
  // But this doesn't have to be the case in future.
  unsigned int nLayers = Utility::countLayers( settings_, stubs_ ); // Count tracker layers with stubs
  bool valid4par = nLayers >= minStubLayers_;
  if (l1track3D.pt() > minPtToReduceLayers_) valid4par = nLayers >= minStubLayers_ - 1;
  bool valid5par = valid4par;
  
  //  Create fitted track objects, specifying track candidate used in the fitting, fitted variables, num of helix parameters used, and Phi sectors and Eta regions.

  if (nPar_ == 4 && valid4par ) return L1fittedTrack(settings_, l1track3D, stubs_, (rinvfit4par_)/(invPtToInvR_), 0., phi0fit4par_, z0fit4par_, tfit4par_, chisq4par_, 4, iPhiSec, iEtaReg, true);
  else if (nPar_ == 5 && valid5par ) return L1fittedTrack(settings_, l1track3D, stubs_, (rinvfit_)/(invPtToInvR_), d0fit_, phi0fit_, z0fit_, tfit_, chisq_, 5, iPhiSec, iEtaReg, true);
  else return L1fittedTrack (settings_, l1track3D, stubs_, l1track3D.qOverPt(), 0., l1track3D.phi0(), l1track3D.z0(), l1track3D.tanLambda(), 999999., 4, iPhiSec, iEtaReg, false);
}

void TrackFitLinearAlgo::invert( float M[5][10],unsigned int n ){
 
  assert(n<=5);
 
  unsigned int i,j,k;
  double ratio,a;
 
  for(i = 0; i < n; i++){
    for(j = n; j < 2*n; j++){
      if(i==(j-n))
        M[i][j] = 1.0;
      else
        M[i][j] = 0.0;
    }
  }
 
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      if(i!=j){
        ratio = M[j][i]/M[i][i];
        for(k = 0; k < 2*n; k++){
          M[j][k] -= ratio * M[i][k];
        }
     }
    }
  }
 
  for(i = 0; i < n; i++){
    a = M[i][i];
    for(j = 0; j < 2*n; j++){
      M[i][j] /= a;
    }
  }
}
 
 
void TrackFitLinearAlgo::calculateDerivatives( bool withd0 ){
 
  unsigned int n=stubs_.size();
 
  //    assert(n<=20);
 
  int j=0;
 
  for(unsigned int i=0;i<n;i++) {
   
    double ri=stubs_[i]->r();
    double zi=stubs_[i]->z();
   
    const double sigmax=stubs_[i]->sigmaX();
    const double sigmaz=stubs_[i]->sigmaZ();

    if ( stubs_[i]->barrel() ){
      //here we handle a barrel hit
 
      //first we have the phi position
      D_[0][j]=-0.5*ri*ri/sqrt(1-0.25*ri*ri*rinv_*rinv_)/sigmax;
      D_[1][j]=ri/sigmax;
      D_[2][j]=0.0;
      D_[3][j]=0.0;
      D_[4][j]=1.0/sigmax;
      j++;
      //second the z position
      D_[0][j]=0.0;
      D_[1][j]=0.0;
      D_[2][j]=(2/rinv_)*asin(0.5*ri*rinv_)/sigmaz;
      D_[3][j]=1.0/sigmaz;
      D_[4][j]=0.0;
      j++;
    }
    else {
      //here we handle a disk hit
      //first we have the r position
 
      double r_track=2.0*sin(0.5*rinv_*(zi-z0_)/t_)/rinv_;
      double phi_track=phi0_-0.5*rinv_*(zi-z0_)/t_;
 
      int iphi=stubs_[i]->iphi();
      double phii=stubs_[i]->phi();
       
      //        std::cout << "iphi: " << iphi << std::endl;
 
      // N.B. These represent HALF the width and number of strips of sensor.
      double width = stubs_[i]->width()/2.0;
      double nstrip = stubs_[i]->nstrip()/2.0;
 
      double Deltai=width*(iphi-nstrip)/nstrip;  //A bit of a hack...
 
      //     std::cout << "Deltai: " << Deltai << std::endl;
      //     std::cout << "ri/zi: " << ri << "/" << stubs_[i]->z() << std::endl;
 
      if (stubs_[i]->z()>0.0) Deltai=-Deltai;
      double theta0=asin(Deltai/ri);
 
      //      std::cout << "theta0: " << theta0 << std::endl;
 
      double rmultiplier=-sin(theta0-(phi_track-phii));
      double phimultiplier=r_track*cos(theta0-(phi_track-phii));
 
 
      double drdrinv=-2.0*sin(0.5*rinv_*(zi-z0_)/t_)/(rinv_*rinv_)
                     +(zi-z0_)*cos(0.5*rinv_*(zi-z0_)/t_)/(rinv_*t_);
      double drdphi0=0;
      double drdt=-(zi-z0_)*cos(0.5*rinv_*(zi-z0_)/t_)/(t_*t_);
      double drdz0=-cos(0.5*rinv_*(zi-z0_)/t_)/t_;
 
      double dphidrinv=-0.5*(zi-z0_)/t_;
      double dphidphi0=1.0;
      double dphidt=0.5*rinv_*(zi-z0_)/(t_*t_);
      double dphidz0=0.5*rinv_/t_;
       
      D_[0][j]=drdrinv/sigmaz;
      D_[1][j]=drdphi0/sigmaz;
      D_[2][j]=drdt/sigmaz;
      D_[3][j]=drdz0/sigmaz;
      D_[4][j]=0.0;
      j++;
      //second the rphi position
      D_[0][j]=(phimultiplier*dphidrinv+rmultiplier*drdrinv)/sigmax;
      D_[1][j]=(phimultiplier*dphidphi0+rmultiplier*drdphi0)/sigmax;
      D_[2][j]=(phimultiplier*dphidt+rmultiplier*drdt)/sigmax;
      D_[3][j]=(phimultiplier*dphidz0+rmultiplier*drdz0)/sigmax;
      D_[4][j]=1.0/sigmax;
      //old calculation
      //D_[0][j]=-0.5*(zi-z0_)/(t_*(sigmax/ri));
      //D_[1][j]=1.0/(sigmax/ri);
      //D_[2][j]=-0.5*rinv_*(zi-z0_)/(t_*t_*(sigmax/ri));
      //D_[3][j]=0.5*rinv_/((sigmax/ri)*t_);
      j++;
    }
 
    //cout << "Exact rinv derivative: "<<i<<" "<<D_[0][j-2]<<" "<<D_[0][j-1]<<endl;
    //cout << "Exact phi0 derivative: "<<i<<" "<<D_[1][j-2]<<" "<<D_[1][j-1]<<endl;
    //cout << "Exact t derivative   : "<<i<<" "<<D_[2][j-2]<<" "<<D_[2][j-1]<<endl;
    //cout << "Exact z0 derivative  : "<<i<<" "<<D_[3][j-2]<<" "<<D_[3][j-1]<<endl;
       
       
  }
   
  //cout << "D:"<<endl;
  //for(unsigned int j=0;j<2*n;j++){
  //  cout <<D_[0][j]<<" "<<D_[1][j]<<" "<<D_[2][j]<<" "<<D_[3][j]<<endl;
  //}
 
     
 
    unsigned int npar=4;
    if ( withd0 ) npar++;
 
    for(unsigned int i1=0;i1<npar;i1++){
      for(unsigned int i2=0;i2<npar;i2++){
        M_[i1][i2]=0.0;
        for(unsigned int j=0;j<2*n;j++){
          M_[i1][i2]+=D_[i1][j]*D_[i2][j];       
        }
      }
    }
 
 
  invert(M_,npar);
 
  for(unsigned int j=0;j<2*n;j++) {
    for(unsigned int i1=0;i1<npar;i1++) {
      MinvDt_[i1][j]=0.0;
      for(unsigned int i2=0;i2<npar;i2++) {
        MinvDt_[i1][j]+=M_[i1][i2+npar]*D_[i2][j];
      }
    }
  }
 
}
 
void TrackFitLinearAlgo::residuals( float& largestresid,int& ilargestresid ) {
 
  unsigned int n=stubs_.size();
 
  //Next calculate the residuals
 
  double delta[2*__MAX_STUBS_PER_TRK__];
 
  double chisq=0.0;
 
  unsigned int j=0;
  
  //  if (print) std::cout << "Residuals ("<<chisq_<<") ["<<invPtToInvR_/rinvfit4par_<<"]: ";
 
  largestresid=-1.0;
  ilargestresid=-1;
 
  for(unsigned int i=0;i<n;i++) {
    double ri=stubs_[i]->r();
    double zi=stubs_[i]->z();
    double phii=stubs_[i]->phi();
    const double sigmax=stubs_[i]->sigmaX();
    const double sigmaz=stubs_[i]->sigmaZ();
 
    if ( stubs_[i]->barrel() ) {
      //we are dealing with a barrel stub
 
      double deltaphi=phi0fit4par_-asin(0.5*ri*rinvfit4par_)-phii;
      if (deltaphi>M_PI) deltaphi-=2*M_PI;
      if (deltaphi<-M_PI) deltaphi+=2*M_PI;
      //      std::cout << "phi0/phii/ri/rinv/ : " << phi0_ << "/" << phii << "/" << ri << "/" << rinv_ << "." << std::endl;
      //      std::cout << "deltaphi: " << deltaphi << std::endl;
      //      assert(fabs(deltaphi)<0.2*M_PI);
 
      delta[j++]=ri*deltaphi/sigmax;
      delta[j++]=(z0fit4par_+(2.0/rinvfit4par_)*tfit4par_*asin(0.5*ri*rinvfit4par_)-zi)/sigmaz;
    }
    else {
      //we are dealing with a disk hit
 
      double r_track=2.0*sin(0.5*rinvfit4par_*(zi-z0fit4par_)/tfit4par_)/rinvfit4par_;
      double phi_track=phi0fit4par_-0.5*rinvfit4par_*(zi-z0fit4par_)/tfit4par_;
       
      int iphi=stubs_[i]->iphi();
 
      // N.B. These represent HALF the width and number of strips of sensor.
      double width = stubs_[i]->width()/2.0;
      double nstrip = stubs_[i]->nstrip()/2.0;
 
      double Deltai=width*(iphi-nstrip)/nstrip;  //A bit of a hack...
 
      if (stubs_[i]->z()>0.0) Deltai=-Deltai;
 
      double theta0=asin(Deltai/ri);
 
      double Delta=Deltai-r_track*sin(theta0-(phi_track-phii));
 
      delta[j++]=(r_track-ri)/sigmaz;
      delta[j++]=Delta/sigmax;
    }
 
    if (fabs(delta[j-2])>largestresid) {
      largestresid=fabs(delta[j-2]);
      ilargestresid=i;
    }
 
    if (fabs(delta[j-1])>largestresid) {
      largestresid=fabs(delta[j-1]);
      ilargestresid=i;
    }
    
    chisq+=delta[j-2]*delta[j-2]+delta[j-1]*delta[j-1];
    if (print) std::cout << __LINE__ << " - Residuals(): delta["<<j-2<<"]/delta["<<j-1<<"]: "<< delta[j-2] << "/" << delta[j-1] << std::endl; 
    if (print) std::cout << __LINE__ << " - Residuals(): chisq: " << chisq  << std::endl; 
 
  }
 
}
 
 
void TrackFitLinearAlgo::linearTrackFit( bool withd0 ) {
 
  unsigned int n=stubs_.size();
 
  //Next calculate the residuals
 
  double delta[2*__MAX_STUBS_PER_TRK__];
 
  double chisq=0;
 
  unsigned int j=0;
 
  for(unsigned int i=0;i<n;i++) {
    double ri=stubs_[i]->r();
    double zi=stubs_[i]->z();
    double phii=stubs_[i]->phi();
    const double sigmax=stubs_[i]->sigmaX();
    const double sigmaz=stubs_[i]->sigmaZ();
 
    if ( stubs_[i]->barrel() ) {
      //we are dealing with a barrel stub
 
      double deltaphi=phi0_-asin(0.5*ri*rinv_)-phii;
      if (deltaphi>M_PI) deltaphi-=2*M_PI;
      if (deltaphi<-M_PI) deltaphi+=2*M_PI;
      //      std::cout << "phi0/phii/ri/rinv/ : " << phi0_ << "/" << phii << "/" << ri << "/" << rinv_ << "." << std::endl;
      //      std::cout << "deltaphi: " << deltaphi << std::endl;
      //      assert(fabs(deltaphi)<0.2*M_PI);
 
      delta[j++]=ri*deltaphi/sigmax;
      delta[j++]=(z0_+(2.0/rinv_)*t_*asin(0.5*ri*rinv_)-zi)/sigmaz;
    }
    else {
      //we are dealing with a disk hit
 
      double r_track=2.0*sin(0.5*rinv_*(zi-z0_)/t_)/rinv_;
      //cout <<"t_track 1: "<<r_track<<endl;
      double phi_track=phi0_-0.5*rinv_*(zi-z0_)/t_;
       
      int iphi=stubs_[i]->iphi();
 
      // N.B. These represent HALF the width and number of strips of sensor.
      double width = stubs_[i]->width()/2.0;
      double nstrip = stubs_[i]->nstrip()/2.0;
 
      double Deltai=width*(iphi-nstrip)/nstrip;  //A bit of a hack...
      //      double Deltai=(stubs_[2].r()-stubs_[1].r()/(1+0.035+3/2*0.035*0.035);  //A bit of a hack...taken to silly levels NB. Not required anymore
       
      if (stubs_[i]->z()>0.0) Deltai=-Deltai;
 
      double theta0=asin(Deltai/ri);
 
      double Delta=Deltai-r_track*sin(theta0-(phi_track-phii));
 
      delta[j++]=(r_track-ri)/sigmaz;
      //double deltaphi=phi_track-phii;
      //if (deltaphi>M_PI) deltaphi-=2*M_PI;
      //if (deltaphi<-M_PI) deltaphi+=2*M_PI;
      //assert(fabs(deltaphi)<0.2*M_PI)
      //delta[j++]=deltaphi/(sigmax/ri);
      delta[j++]=Delta/sigmax;
 
    }
 
    chisq+=(delta[j-2]*delta[j-2]+delta[j-1]*delta[j-1]);
    if (print) std::cout << __LINE__ << " - linearTrackFit(): delta["<<j-2<<"]/delta["<<j-1<<"] "<< delta[j-2] << "/" << delta[j-1] << std::endl; 
  }
 
  double drinv = 0.0;
  double dphi0 = 0.0;
  double dd0 = 0.0;
  double dt = 0.0;
  double dz0 = 0.0;
 
  double drinv_cov = 0.0;
  double dphi0_cov = 0.0;
  double dd0_cov = 0.0;
  double dt_cov = 0.0;
  double dz0_cov = 0.0;
 
 
 
  for(unsigned int j=0;j<2*n;j++) {
    drinv-=MinvDt_[0][j]*delta[j];
    dphi0-=MinvDt_[1][j]*delta[j];
    dt-=MinvDt_[2][j]*delta[j];
    dz0-=MinvDt_[3][j]*delta[j];
    if ( withd0 ) dd0-=MinvDt_[4][j]*delta[j];
 
    drinv_cov+=D_[0][j]*delta[j];
    dphi0_cov+=D_[1][j]*delta[j];
    dt_cov+=D_[2][j]*delta[j];
    dz0_cov+=D_[3][j]*delta[j];
    if ( withd0 ) dd0_cov+=D_[4][j]*delta[j];

  }
   
  //std::cout << __LINE__ << " - linearTrackFit(): chisq: "<< chisq << std::endl;

  double deltaChisq=drinv*drinv_cov+dphi0*dphi0_cov+dt*dt_cov+dz0*dz0_cov;
  if ( withd0 ) deltaChisq+=dd0*dd0_cov;
  //std::cout << __LINE__ << "- deltaChisq: " << deltaChisq << std::endl; 
  //std::cout << __LINE__ << " - linearTrackFit(): chisq+deltaChisq: " << chisq+deltaChisq << std::endl;

  //drinv=0.0; dphi0=0.0; dt=0.0; dz0=0.0;
 
    if ( withd0 ) {
      rinvfit_=rinv_+drinv;
      phi0fit_=phi0_+dphi0;
     
      tfit_=t_+dt;
      z0fit_=z0_+dz0;
     
      d0fit_=d0_+dd0;
 
      chisq_=(chisq+deltaChisq);
 
    }
   
    else {
      rinvfit4par_=rinv_+drinv;
      phi0fit4par_=phi0_+dphi0;
      tfit4par_=t_+dt;
      z0fit4par_=z0_+dz0;
     
      chisq4par_=(chisq+deltaChisq);
    }
 
  //cout << "Trackfit:"<<endl;
  //cout << "rinv_ drinv: "<<rinv_<<" "<<drinv<<endl;
  //cout << "phi0_ dphi0: "<<phi0_<<" "<<dphi0<<endl;
  //cout << "t_ dt      : "<<t_<<" "<<dt<<endl;
  //cout << "z0_ dz0    : "<<z0_<<" "<<dz0<<endl;
 
}
