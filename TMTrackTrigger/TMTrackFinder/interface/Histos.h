#ifndef __HISTOS_H__
#define __HISTOS_H__

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TMTrackTrigger/TMTrackFinder/interface/Settings.h>

#include "boost/numeric/ublas/matrix.hpp"

#include <vector>
#include <map>
#include <string>


// #define HISTOS_OPTIMIZE_


class InputData;
class TP;
class Sector;
class HTpair;
class L1fittedTrack;
class L1fittedTrk4and5;
class TH1F;
class TH2F;
class TProfile;
class TGraphAsymmErrors;

class Histos {

public:
	// Store cfg parameters.
	Histos(const Settings* settings) : settings_(settings), numPerfRecoTPforAlg_(0) {}

	~Histos(){}

	// Book all histograms
	void book();

	// Fill histograms with stubs and tracking particles from input data.
	void fillInputData(const InputData& inputData);
	// Fill histograms that check if choice of (eta,phi) sectors is good.
	void fillEtaPhiSectors(const InputData& inputData, const boost::numeric::ublas::matrix<Sector>& mSectors);
	// Fill histograms checking filling of r-phi HT array.
	void fillRphiHT(const boost::numeric::ublas::matrix<HTpair>& mHtPairs);
	// Book histograms about r-z track filters (or other filters applied after r-phi HT array).
	void fillRZfilters(const boost::numeric::ublas::matrix<HTpair>& mHtPairs);
	// Fill histograms studying track candidates found by Hough Transform.
	void fillTrackCands(const InputData& inputData, const boost::numeric::ublas::matrix<Sector>& mSectors, const boost::numeric::ublas::matrix<HTpair>& mHtPairs);
	// Fill histograms studying freak, events with too many stubs..
	void fillStudyBusyEvents(const InputData& inputData, const boost::numeric::ublas::matrix<Sector>& mSectors, const boost::numeric::ublas::matrix<HTpair>& mHtPairs);
	// Fill histograms relating to track fitting performance.
	void fillTrackFitting(const InputData& inputData, const std::vector<std::pair<std::string,L1fittedTrack>>& fittedTracks, float chi2dofCutPlots);

	void endJobAnalysis();

private:

	// Book histograms for specific topics.
	void bookInputData();
	void bookEtaPhiSectors();
	void bookRphiHT();
	void bookRZfilters();
	void bookTrackCands();
	void bookStudyBusyEvents();
	void bookTrackFitting();

	// Produce plots of tracking efficiency prior to track fit (run at end of job).
	void plotTrackEfficiency();
	// Produce plots of tracking efficiency after track fit (run at end of job).
	void plotTrackEffAfterFit(std::string fitName);

	// Understand why not all tracking particles were reconstructed.
	// Returns list of tracking particles that were not reconstructed and an integer indicating why.
	// Only considers TP used for algorithmic efficiency measurement.
	std::map<const TP*, std::string> diagnoseTracking(const InputData& inputData, const boost::numeric::ublas::matrix<Sector>& mSectors, const boost::numeric::ublas::matrix<HTpair>& mHtPairs) const;

private:

	const Settings *settings_; // Configuration parameters.

	edm::Service<TFileService> fs_;

	// Histograms of input data.
	TProfile* profNumStubs_;
	TH1F* hisStubsVsEta_;
	TH1F* hisStubsVsR_;
	TProfile* profNumTPs_;
	TH1F* hisNumStubsPerTP_;
	TProfile* hisStubKillFE_;
	TProfile* hisStubIneffiVsInvPt_;
	TProfile* hisStubIneffiVsEta_;
	TProfile* hisStubKillDataCorr_;
	TH1F* hisPtStub_;
	TH1F* hisPtResStub_;
	TH1F* hisBendFilterPower_;
	TH1F* hisDelPhiStub_;
	TH1F* hisDelPhiResStub_;
	TH1F* hisBendStub_;
	TH1F* hisBendResStub_;
	TH1F* hisNumMergedBend_;
	TH2F* hisBendVsLayerOrRing_;
	TH2F* hisBendFEVsLayerOrRing_;
	TH1F* hisPhiStubVsPhiTP_;
	TH1F* hisPhiStubVsPhi0TP_;
	TH1F* hisPhi0StubVsPhi0TP_;
	TH1F* hisPhi0StubVsPhi0TPres_;
	TH1F* hisPhiStubVsPhi65TP_;
	TH1F* hisPhi65StubVsPhi65TP_;
	TH1F* hisPhi65StubVsPhi65TPres_;
	TH1F* hisPitchOverSep_;
	TH1F* hisRhoParameter_;
	TH1F* hisFracStubsSharingClus0_;
	TH1F* hisFracStubsSharingClus1_;

	// Histograms checking that (eta,phi) sector definition is good.
	TH1F* hisFracStubsInSec_;
	TH1F* hisFracStubsInEtaSec_;
	TH1F* hisFracStubsInPhiSec_;
	TH1F* hisNumSecsPerStub_;
	TH1F* hisNumEtaSecsPerStub_;
	TH1F* hisNumPhiSecsPerStub_;
	TH1F* hisNumStubsPerSec_;
	TProfile* profNumStubsPerEtaSec_;
	TH2F* hisLayerIDvsEtaSec_;
	TH2F* hisLayerIDreducedvsEtaSec_;

	// Histograms checking filling of r-phi HT array.
	TH1F* hisIncStubsPerHT_;
	TH1F* hisExcStubsPerHT_;
	TH2F* hisNumStubsInCellVsEta_;
	TH1F* hisStubsOnRphiTracksPerHT_;

	// Histograms about r-z track filters (or other filters applied after r-phi HT array).
	TH1F* hisNumZtrkSeedCombinations_;
	TH1F* hisNumSeedCombinations_;
	TH1F* hisNumGoodSeedCombinations_;
	TH1F* hisCorrelationZTrk_;

	// Histograms studying track candidates found by Hough Transform.
	TProfile* profNumTrackCands_;
	TProfile* profNumTracksVsEta_;
	TH1F*     hisNumTracksVsQoverPt_;
	TH1F*     hisNumTrksPerSect_;
	TH1F*     hisNumTrksPerOct_;
	TProfile* profStubsOnTracks_;
	TProfile* profStubsOnTracksVsEta_;
	TH1F*     hisStubsOnTracksPerSect_;
	TH1F*     hisStubsOnTracksPerOct_;
	TH1F*     hisUniqueStubsOnTrksPerSect_;
	TH1F*     hisUniqueStubsOnTrksPerOct_;
	TH1F*     hisStubsPerTrack_;
	TH1F*     hisLayersPerTrack_;
	TH1F*     hisPSLayersPerTrack_;
	TH1F*     hisLayersPerTrueTrack_;
	TH1F*     hisPSLayersPerTrueTrack_;
	TProfile* profExcessStubsPerTrackVsPt_;
	TH1F*     hisFracMatchStubsOnTracks_;
	TH1F* hisDeltaPhiRtruePS_;
	TH1F* hisDeltaRorZtruePS_;
	TH1F* hisDeltaPhiRtrue2S_;
	TH1F* hisDeltaRorZtrue2S_;
	TH1F* hisDeltaPhiRfakePS_;
	TH1F* hisDeltaRorZfakePS_;
	TH1F* hisDeltaPhiRfake2S_;
	TH1F* hisDeltaRorZfake2S_;
	TProfile* profNsigmaPhiRvsInvPt_;
	TProfile* profNsigmaPhiRvsFracDist_;
	TH1F* hisDeltaBendTrue_;
	TH1F* hisDeltaBendFake_;
	TProfile* profFracTrueStubsVsLayer_;
	TProfile* profDupTracksVsTPeta_;

	// Histograms of track parameter resolution after HT transform.
	TH1F* hisQoverPtRes_;
	TH1F* hisPhi0Res_;
	TH1F* hisEtaRes_;
	TH1F* hisZ0Res_;

	// Diagnosis of failed tracking.
	TH1F* hisRecoFailureReason_;
	TH1F* hisRecoFailureLayer_;

	// Histograms used to make efficiency plots with track candidates prior to fit.
	TH1F* hisTPinvptForEff_;
	TH1F* hisRecoTPinvptForEff_;
	TH1F* hisTPetaForEff_;
	TH1F* hisRecoTPetaForEff_;
	TH1F* hisTPphiForEff_;
	TH1F* hisRecoTPphiForEff_;
	//
	TH1F* hisPerfRecoTPinvptForEff_;
	//
	TH1F* hisTPinvptForAlgEff_;
	TH1F* hisRecoTPinvptForAlgEff_;
	TH1F* hisTPetaForAlgEff_;
	TH1F* hisRecoTPetaForAlgEff_;
	TH1F* hisTPphiForAlgEff_;
	TH1F* hisRecoTPphiForAlgEff_;
	//
	TH1F* hisPerfRecoTPinvptForAlgEff_;
	//
	TH1F* hisTPd0ForAlgEff_;
	TH1F* hisRecoTPd0ForAlgEff_;
	TH1F* hisTPz0ForAlgEff_;
	TH1F* hisRecoTPz0ForAlgEff_;

	// Histograms for studying freak, large events with too many stubs.
	TH1F*     hisNumBusySecsInPerEvent_;
	TH1F*     hisNumBusySecsOutPerEvent_;
	TProfile* profFracBusyInVsEtaReg_;
	TProfile* profFracBusyOutVsEtaReg_;
	TProfile* profFracStubsKilledVsEtaReg_;
	TProfile* profFracTracksKilledVsEtaReg_;
	TProfile* profFracTracksKilledVsInvPt_;
	TProfile* profFracTPKilledVsEta_;
	TProfile* profFracTPKilledVsInvPt_;
	TH1F*     hisNumTPkilledBusySec_;
	
	
#ifndef HISTOS_OPTIMIZE_
	
	std::map<std::string, TH1F*> hisNumInputStubs_;
	std::map<std::string, TH1F*> hisQoverPtInputStubs_;
	std::map<std::string, TH1F*> hisNumOutputStubs_;
	std::map<std::string, TH1F*> hisNumTracks_; 
	std::map<std::string, TH1F*> hisNumStubsPerTrack_; 
	std::map<std::string, TH1F*> hisTrackQoverPt_; 
	std::map<std::string, TH1F*> hisTrackPurity_; 
	std::map<std::string, TH1F*> hisNumTPphysics_; 
	std::map<std::string, TH1F*> hisNumTPpileup_; 
	std::map<std::string, TH1F*> hisSumPtTPphysics_; 
	std::map<std::string, TH1F*> hisSumPtTPpileup_; 
	
#else
	
	enum ttype {BusyOutSec, QuietOutSec, NumberOfTNames};
	std::array<TH1F*, NumberOfTNames> hisNumInputStubs_;
	std::array<TH1F*, NumberOfTNames> hisQoverPtInputStubs_;
	std::array<TH1F*, NumberOfTNames> hisNumOutputStubs_;
	std::array<TH1F*, NumberOfTNames> hisNumTracks_; 
	std::array<TH1F*, NumberOfTNames> hisNumStubsPerTrack_; 
	std::array<TH1F*, NumberOfTNames> hisTrackQoverPt_; 
	std::array<TH1F*, NumberOfTNames> hisTrackPurity_; 
	std::array<TH1F*, NumberOfTNames> hisNumTPphysics_; 
	std::array<TH1F*, NumberOfTNames> hisNumTPpileup_; 
	std::array<TH1F*, NumberOfTNames> hisSumPtTPphysics_; 
	std::array<TH1F*, NumberOfTNames> hisSumPtTPpileup_;
	
#endif


#ifndef HISTOS_OPTIMIZE_

	// Histograms for track fitting evaluation, where map index specifies name of track fitting algorithm used.
	std::map<std::string, TH1F*    > hisSeedQinvPt_;
	std::map<std::string, TH1F*    > hisSeedPhi0_;
	std::map<std::string, TH1F*    > hisSeedD0_;
	std::map<std::string, TH1F*    > hisSeedZ0_;
	std::map<std::string, TH1F*    > hisSeedEta_;

	std::map<std::string, TProfile*> profNumFittedCands_;

	std::map<std::string, TH1F*    > hisFitQinvPtMatched_;
	std::map<std::string, TH1F*    > hisFitPhi0Matched_;
	std::map<std::string, TH1F*    > hisFitD0Matched_;
	std::map<std::string, TH1F*    > hisFitZ0Matched_;
	std::map<std::string, TH1F*    > hisFitEtaMatched_;

	std::map<std::string, TH1F*    > hisFitChi2Matched_;
	std::map<std::string, TH1F*    > hisFitChi2DofMatched_;

	std::map<std::string, TH1F*    > hisFitQinvPtUnmatched_;
	std::map<std::string, TH1F*    > hisFitPhi0Unmatched_;
	std::map<std::string, TH1F*    > hisFitD0Unmatched_;
	std::map<std::string, TH1F*    > hisFitZ0Unmatched_;
	std::map<std::string, TH1F*    > hisFitEtaUnmatched_;

	std::map<std::string, TH1F*    > hisFitChi2Unmatched_;
	std::map<std::string, TH1F*    > hisFitChi2DofUnmatched_;

	std::map<std::string, TH2F*    > hisFitVsTrueQinvPtGoodChi2_;
	std::map<std::string, TH2F*    > hisFitVsTruePhi0GoodChi2_;
	std::map<std::string, TH2F*    > hisFitVsTrueD0GoodChi2_;
	std::map<std::string, TH2F*    > hisFitVsTrueZ0GoodChi2_;
	std::map<std::string, TH2F*    > hisFitVsTrueEtaGoodChi2_;

	std::map<std::string, TH2F*    > hisFitVsSeedQinvPtGenCand_;
	std::map<std::string, TH2F*    > hisFitVsSeedPhi0GenCand_;
	std::map<std::string, TH2F*    > hisFitVsSeedD0GenCand_;
	std::map<std::string, TH2F*    > hisFitVsSeedZ0GenCand_;
	std::map<std::string, TH2F*    > hisFitVsSeedEtaGenCand_;

	std::map<std::string, TH1F*    > hisFitQinvPtResGoodChi2_;
	std::map<std::string, TH1F*    > hisFitPhi0ResGoodChi2_;
	std::map<std::string, TH1F*    > hisFitD0ResGoodChi2_;
	std::map<std::string, TH1F*    > hisFitZ0ResGoodChi2_;
	std::map<std::string, TH1F*    > hisFitEtaResGoodChi2_;  

	std::map<std::string, TH1F*    > hisSeedQinvPtResGoodChi2_;
	std::map<std::string, TH1F*    > hisSeedPhi0ResGoodChi2_;
	std::map<std::string, TH1F*    > hisSeedD0ResGoodChi2_;
	std::map<std::string, TH1F*    > hisSeedZ0ResGoodChi2_;
	std::map<std::string, TH1F*    > hisSeedEtaResGoodChi2_;  

	std::map<std::string, TH2F*    > hisFitVsSeedQinvPtFakeCand_;
	std::map<std::string, TH2F*    > hisFitVsSeedPhi0FakeCand_;
	std::map<std::string, TH2F*    > hisFitVsSeedD0FakeCand_;
	std::map<std::string, TH2F*    > hisFitVsSeedZ0FakeCand_;
	std::map<std::string, TH2F*    > hisFitVsSeedEtaFakeCand_;

	std::map<std::string, TProfile*> hisQoverPtResVsTrueEta_;
	std::map<std::string, TProfile*> hisPhi0ResVsTrueEta_;
	std::map<std::string, TProfile*> hisEtaResVsTrueEta_;
	std::map<std::string, TProfile*> hisZ0ResVsTrueEta_;
	std::map<std::string, TProfile*> hisD0ResVsTrueEta_;

	std::map<std::string, TProfile*> hisQoverPtResVsTrueInvPt_;
	std::map<std::string, TProfile*> hisPhi0ResVsTrueInvPt_;
	std::map<std::string, TProfile*> hisEtaResVsTrueInvPt_;
	std::map<std::string, TProfile*> hisZ0ResVsTrueInvPt_;
	std::map<std::string, TProfile*> hisD0ResVsTrueInvPt_;

	std::map<std::string, TH2F*    > hisTrueFittedChiSquaredVsTrueEta_;
	std::map<std::string, TH2F*    > hisTrueFittedChiSquaredDofVsTrueEta_;
	std::map<std::string, TH2F*    > hisTrueFittedChiSquaredVsFittedEta_;
	std::map<std::string, TH2F*    > hisTrueFittedChiSquaredDofVsFittedEta_;
	std::map<std::string, TH2F*    > hisFittedChiSquaredFunctionOfStubs_;
	std::map<std::string, TH2F*    > hisFittedChiSquaredDofFunctionOfStubs_;

	std::map<std::string, TH1F*    > hisTrueEtaMatchedGoodChi2_;
	std::map<std::string, TH1F*    > hisTrueEtaMatchedBadChi2_;
	std::map<std::string, TH1F*    > hisStubPurityMatchedGoodChi2_;
	std::map<std::string, TH1F*    > hisStubPurityMatchedBadChi2_;

	std::map<std::string, TProfile*> profChi2DofVsInvPtPERF_;
	std::map<std::string, TProfile*> profBigChi2DofVsInvPtPERF_;
	std::map<std::string, TH1F*    > hisD0TPBigChi2DofPERF_;
	std::map<std::string, TH1F*    > hisD0TPSmallChi2DofPERF_;

	std::map<std::string, TH2F*    > hisNumMatchedStubsKilledVsKilled_;
	std::map<std::string, TProfile*> profTrksKilledByFit_;
	std::map<std::string, TH2F*    > hisNumStubsVsPurity_;

	std::map<std::string, TH1F*    > hisNumFittingIterations_;
	std::map<std::string, TH2F*    > hisNumFittingIterationsVsPurity_;
	std::map<std::string, TH2F*    > hisNumFittingIterationsVsPurityMatched_;
	std::map<std::string, TH2F*    > hisNumFittingIterationsVsPurityUnmatched_;

	std::map<std::string, TH2F*    > hisFitEfficiencyVsChi2Dof_;
	std::map<std::string, TH2F*    > hisNumStubsVsChi2Dof_;
	std::map<std::string, TH2F*    > hisNumLayersVsChi2Dof_;
	std::map<std::string, TH2F*    > hisAvgNumStubsPerLayerVsChi2Dof_;

	// Histograms used for efficiency plots made with fitted tracks.
	std::map<std::string, TH1F*    > hisFitTPinvptForEff_;
	std::map<std::string, TH1F*    > hisFitTPetaForEff_;
	std::map<std::string, TH1F*    > hisFitTPphiForEff_;
	std::map<std::string, TH1F*    > hisPerfFitTPinvptForEff_;
	std::map<std::string, TH1F*    > hisFitTPinvptForAlgEff_;
	std::map<std::string, TH1F*    > hisFitTPetaForAlgEff_;
	std::map<std::string, TH1F*    > hisFitTPphiForAlgEff_;
	std::map<std::string, TH1F*    > hisPerfFitTPinvptForAlgEff_;
	std::map<std::string, TH1F*    > hisFitTPd0ForAlgEff_;
	std::map<std::string, TH1F*    > hisFitTPz0ForAlgEff_;

#else
	
	enum {MAX_NUMBER_OF_FITTERS = 8};
	std::unordered_map<std::string, std::size_t> fitterNameToFitterIndexMap_;
	
	// Histograms for track fitting evaluation, where map index specifies name of track fitting algorithm used.
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisSeedQinvPt_                            ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisSeedPhi0_                              ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisSeedD0_                                ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisSeedZ0_                                ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisSeedEta_                               ;

	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> profNumFittedCands_                       ;

	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitQinvPtMatched_                      ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitPhi0Matched_                        ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitD0Matched_                          ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitZ0Matched_                          ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitEtaMatched_                         ;

	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitChi2Matched_                        ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitChi2DofMatched_                     ;

	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitQinvPtUnmatched_                    ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitPhi0Unmatched_                      ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitD0Unmatched_                        ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitZ0Unmatched_                        ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitEtaUnmatched_                       ;

	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitChi2Unmatched_                      ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitChi2DofUnmatched_                   ;

	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsTrueQinvPtGoodChi2_               ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsTruePhi0GoodChi2_                 ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsTrueD0GoodChi2_                   ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsTrueZ0GoodChi2_                   ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsTrueEtaGoodChi2_                  ;

	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsSeedQinvPtGenCand_                ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsSeedPhi0GenCand_                  ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsSeedD0GenCand_                    ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsSeedZ0GenCand_                    ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsSeedEtaGenCand_                   ;

	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitQinvPtResGoodChi2_                  ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitPhi0ResGoodChi2_                    ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitD0ResGoodChi2_                      ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitZ0ResGoodChi2_                      ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitEtaResGoodChi2_                     ;

	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisSeedQinvPtResGoodChi2_                 ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisSeedPhi0ResGoodChi2_                   ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisSeedD0ResGoodChi2_                     ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisSeedZ0ResGoodChi2_                     ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisSeedEtaResGoodChi2_                    ;

	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsSeedQinvPtFakeCand_               ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsSeedPhi0FakeCand_                 ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsSeedD0FakeCand_                   ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsSeedZ0FakeCand_                   ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitVsSeedEtaFakeCand_                  ;

	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> hisQoverPtResVsTrueEta_                   ;
	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> hisPhi0ResVsTrueEta_                      ;
	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> hisEtaResVsTrueEta_                       ;
	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> hisZ0ResVsTrueEta_                        ;
	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> hisD0ResVsTrueEta_                        ;

	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> hisQoverPtResVsTrueInvPt_                 ;
	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> hisPhi0ResVsTrueInvPt_                    ;
	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> hisEtaResVsTrueInvPt_                     ;
	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> hisZ0ResVsTrueInvPt_                      ;
	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> hisD0ResVsTrueInvPt_                      ;

	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisTrueFittedChiSquaredVsTrueEta_         ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisTrueFittedChiSquaredDofVsTrueEta_      ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisTrueFittedChiSquaredVsFittedEta_       ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisTrueFittedChiSquaredDofVsFittedEta_    ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFittedChiSquaredFunctionOfStubs_       ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFittedChiSquaredDofFunctionOfStubs_    ;

	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisTrueEtaMatchedGoodChi2_                ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisTrueEtaMatchedBadChi2_                 ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisStubPurityMatchedGoodChi2_             ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisStubPurityMatchedBadChi2_              ;

	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> profChi2DofVsInvPtPERF_                   ;
	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> profBigChi2DofVsInvPtPERF_                ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisD0TPBigChi2DofPERF_                    ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisD0TPSmallChi2DofPERF_                  ;

	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisNumMatchedStubsKilledVsKilled_         ;
	std::array<TProfile*  , MAX_NUMBER_OF_FITTERS> profTrksKilledByFit_                      ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisNumStubsVsPurity_                      ;

	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisNumFittingIterations_                  ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisNumFittingIterationsVsPurity_          ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisNumFittingIterationsVsPurityMatched_   ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisNumFittingIterationsVsPurityUnmatched_ ;

	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisFitEfficiencyVsChi2Dof_                ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisNumStubsVsChi2Dof_                     ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisNumLayersVsChi2Dof_                    ;
	std::array<TH2F*      , MAX_NUMBER_OF_FITTERS> hisAvgNumStubsPerLayerVsChi2Dof_          ;

	// Histograms used for efficiency plots made with fitted tracks.
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitTPinvptForEff_                      ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitTPetaForEff_                        ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitTPphiForEff_                        ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisPerfFitTPinvptForEff_                  ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitTPinvptForAlgEff_                   ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitTPetaForAlgEff_                     ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitTPphiForAlgEff_                     ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisPerfFitTPinvptForAlgEff_               ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitTPd0ForAlgEff_                      ;
	std::array<TH1F*      , MAX_NUMBER_OF_FITTERS> hisFitTPz0ForAlgEff_                      ;
	
#endif



	// Histograms of tracking efficiency & fake rate after Hough transform based on tracks prior to track fit.
	TGraphAsymmErrors* graphEffVsInvPt_;
	TGraphAsymmErrors* graphEffVsEta_;
	TGraphAsymmErrors* graphEffVsPhi_;
	//
	TGraphAsymmErrors* graphPerfEffVsInvPt_;
	//
	TGraphAsymmErrors* graphAlgEffVsInvPt_;
	TGraphAsymmErrors* graphAlgEffVsEta_;
	TGraphAsymmErrors* graphAlgEffVsPhi_;
	//
	TGraphAsymmErrors* graphPerfAlgEffVsInvPt_;
	//
	TGraphAsymmErrors* graphAlgEffVsD0_;
	TGraphAsymmErrors* graphAlgEffVsZ0_;
	

	// Number of perfectly reconstructed tracks amongst TP used for algorithmic efficiency measurement.
	// Perfectly means that all stubs on track were produced by same TP.
	unsigned int numPerfRecoTPforAlg_;
	
#ifndef HISTOS_OPTIMIZE_
	
	// Histograms of tracking efficiency & fake rate after Hough transform based on tracks after the track fit.
	std::map<std::string, TGraphAsymmErrors*> graphEffFitVsInvPt_;
	std::map<std::string, TGraphAsymmErrors*> graphEffFitVsEta_;
	std::map<std::string, TGraphAsymmErrors*> graphEffFitVsPhi_;
	//
	std::map<std::string, TGraphAsymmErrors*> graphPerfEffFitVsInvPt_;
	//
	std::map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsInvPt_;
	std::map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsEta_;
	std::map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsPhi_;
	//
	std::map<std::string, TGraphAsymmErrors*> graphPerfAlgEffFitVsInvPt_;
	//
	std::map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsD0_;
	std::map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsZ0_;
	

	// Number of genuine reconstructed and perfectly reconstructed tracks which were fitted.
	std::map<std::string, unsigned int> numFitAlgEff_;
	std::map<std::string, unsigned int> numFitPerfAlgEff_;

	// Number of genuine reconstructed and perfectly reconstructed tracks which were fitted post-cut.
	std::map<std::string, unsigned int> numFitAlgEffPass_;
	std::map<std::string, unsigned int> numFitPerfAlgEffPass_;
	
#else
	
	// Histograms of tracking efficiency & fake rate after Hough transform based on tracks after the track fit.
	std::array<TGraphAsymmErrors*, MAX_NUMBER_OF_FITTERS> graphEffFitVsInvPt_;
	std::array<TGraphAsymmErrors*, MAX_NUMBER_OF_FITTERS> graphEffFitVsEta_;
	std::array<TGraphAsymmErrors*, MAX_NUMBER_OF_FITTERS> graphEffFitVsPhi_;
	//
	std::array<TGraphAsymmErrors*, MAX_NUMBER_OF_FITTERS> graphPerfEffFitVsInvPt_;
	//
	std::array<TGraphAsymmErrors*, MAX_NUMBER_OF_FITTERS> graphAlgEffFitVsInvPt_;
	std::array<TGraphAsymmErrors*, MAX_NUMBER_OF_FITTERS> graphAlgEffFitVsEta_;
	std::array<TGraphAsymmErrors*, MAX_NUMBER_OF_FITTERS> graphAlgEffFitVsPhi_;
	//
	std::array<TGraphAsymmErrors*, MAX_NUMBER_OF_FITTERS> graphPerfAlgEffFitVsInvPt_;
	//
	std::array<TGraphAsymmErrors*, MAX_NUMBER_OF_FITTERS> graphAlgEffFitVsD0_;
	std::array<TGraphAsymmErrors*, MAX_NUMBER_OF_FITTERS> graphAlgEffFitVsZ0_;
	

	// Number of genuine reconstructed and perfectly reconstructed tracks which were fitted.
	std::array<unsigned int, MAX_NUMBER_OF_FITTERS> numFitAlgEff_;
	std::array<unsigned int, MAX_NUMBER_OF_FITTERS> numFitPerfAlgEff_;

	// Number of genuine reconstructed and perfectly reconstructed tracks which were fitted post-cut.
	std::array<unsigned int, MAX_NUMBER_OF_FITTERS> numFitAlgEffPass_;
	std::array<unsigned int, MAX_NUMBER_OF_FITTERS> numFitPerfAlgEffPass_;
	
#endif
	
};

#endif
