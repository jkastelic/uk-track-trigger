#include "FWCore/Utilities/interface/Exception.h"

#include "TMTrackTrigger/TMTrackFinder/interface/TrkRZfilter.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"

//=== Initialize configuration parameters, and note eta range covered by sector and phi coordinate of its centre.

void TrkRZfilter::init(const Settings* settings, float etaMinSector, float etaMaxSector, float phiCentreSector) {

  // Configuration parameters.
  settings_ = settings;

  // Eta range of sector & phi coord of its centre.
  etaMinSector_ = etaMinSector;
  etaMaxSector_ = etaMaxSector;
  phiCentreSector_ = phiCentreSector;

  // Calculate z coordinate a track would have at given radius if on eta sector boundary.
  chosenRofZ_    = settings->chosenRofZ();
  zTrkMinSector_ = chosenRofZ_/ tan( 2. * atan(exp(-etaMinSector_)) );
  zTrkMaxSector_ = chosenRofZ_/ tan( 2. * atan(exp(-etaMaxSector_)) );

  // Use filter in each HT cell using only stubs which have consistent rapidity?
  useEtaFilter_  = settings->useEtaFilter();
  // Use filter in each HT cell using only stubs which have consistent zR
  useZTrkFilter_ = settings->useZTrkFilter();
  // Filter stubs in cell using a tracklet-like algorithm
  useSeedFilter_ = settings->useSeedFilter();

  // --- Options for ZTrk filter
  // Use z of track at this radius for ZTrkFilter.
  chosenRofZFilter_    = settings->chosenRofZFilter();

  // --- Options for Seed filter.
  //Added resolution for a tracklet-like filter algorithm, beyond that estimated from hit resolution.
  seedResolution_      = settings->seedResolution();
  // Keep stubs compatible with all possible good seed.
  keepAllSeed_         = settings->keepAllSeed();
  // Maximum number of seed combinations to bother checking per track candidate.
  maxSeedCombinations_ = settings->maxSeedCombinations();
  // Reject tracks whose estimated rapidity from seed filter is inconsistent range of with eta sector. (Kills some duplicate tracks).
  zTrkSectorCheck_     = settings->zTrkSectorCheck();

  // Min. number of layers track must have stubs in to be declared valid.
  minStubLayers_       = settings->minStubLayers();
  // Reduce MinStubLayers by 1 for tracks with pt above this cut.
  // If this is set to >= 10000., this option is disabled.
  minPtToReduceLayers_ = settings->minPtToReduceLayers();

  // Assumed length of beam-spot in z.
  beamWindowZ_ = settings->beamWindowZ();

  // Note that no r-z filter has yet provided an estimate of the r-z track parameters.
  estValid_ = false;  // No valid estimate yet.
}

// Filters track candidates (found by the r-phi Hough transform), removing inconsistent stubs from the tracks, 
// also killing some of the tracks altogether if they are left with too few stubs.
// Also adds an estimate of r-z helix parameters to the selected track objects, if the filters used provide this.
 
vector<L1track2D> TrkRZfilter::filterTracks(const vector<L1track2D>& tracks) {

  vector<L1track2D> filteredTracks;

  for (const L1track2D& trkIN : tracks) {
    if (! trkIN.isRphiTrk()) throw cms::Exception("TrkRZfilter ERROR: Called for track found by r-z Hough transform!");    

    const vector<const Stub*>& stubs = trkIN.getStubs(); // stubs assigned to track 

    // Filter stubs assigned to track, checking they are consistent with requested criteria.
    vector<const Stub*> filteredStubs = stubs;
    if (useEtaFilter_)  filteredStubs = this->etaFilter (filteredStubs, trkIN.qOverPt());
    if (useZTrkFilter_) filteredStubs = this->zTrkFilter(filteredStubs, trkIN.qOverPt());
    if (useSeedFilter_) filteredStubs = this->seedFilter(filteredStubs, trkIN.qOverPt());

    // Check if track still has stubs in enough layers after filter.
    unsigned int numLayersAfterFilters = Utility::countLayers(settings_, filteredStubs);
    if ( this->trackCandCheck( numLayersAfterFilters, trkIN.qOverPt() ) ) {

      // Create copy of original track, except now using its filtered stubs, to be added to filteredTrack collection.
      L1track2D trkOUT(settings_, filteredStubs, trkIN.getCellLocation(), trkIN.getHelix2D(), true);

      // If one of the filters provided an estimate of the r-z track parameters, then store it inside track.
      if (estValid_) trkOUT.setTrkEstZ0andTanLam(estZ0_, estTanLambda_);

      filteredTracks.push_back(trkOUT);
    }
  }

  return filteredTracks;
}

//=== Produce a filtered collection of stubs on this track candidate that all have consistent rapidity.
//=== Only called for r-phi Hough transform.

// Davide's algorithm
vector<const Stub*> TrkRZfilter::etaFilter( const vector<const Stub*>& stubs, float trkQoverPt ) const {

  // Define histogram in stub rapidity. (Should improve this by using sector rapidity range!).
  const unsigned int nBins = 64;
  const float etaMax =  3.1;
  const float etaMin = -etaMax;
  const float etaBinWidth = (etaMax - etaMin)/nBins;
  vector<unsigned int> eta_hist(nBins, 0);

  // Loop over stubs, histogramming their rapidity.
  for (const Stub* s : stubs) {
    int ibin = std::floor( (s->eta() - etaMin) / etaBinWidth );
    if (ibin >= 0 && ibin < int(nBins)) eta_hist[ibin]++; 
  }

  // Find mode of histogram.
  int modeBin = -1;
  int maxEntry = -1;  
  for (unsigned int i = 0; i < nBins; i++) {
    if (maxEntry < int(eta_hist[i])) {
      modeBin = i;
      maxEntry = eta_hist[i];
    }
  }

  float etaBest = etaMin + (modeBin + 0.5)*etaBinWidth; // Centre of mode bin.
  float etaResolution = (-0.0775 * fabs(etaBest)) + 0.35; // assumed resolution in eta
  //  float etaResolution = 0.5*0.15*sin(2*atan(exp(-2*etaBest)))*2 + 0.05; // assumed resolution in eta + half bin width contribution.

  // Create eta-filtered stub collection.
  vector<const Stub*> filteredStubs;
  for (const Stub* s : stubs) {
    if (fabs(s->eta() - etaBest) < etaResolution) filteredStubs.push_back(s);
  }

  return filteredStubs;
}

//=== Produce a filtered collection of stubs on this track candidate that all have consistent zR.

vector<const Stub*> TrkRZfilter::zTrkFilter(const vector<const Stub*>& stubs, float trkQoverPt ) {
  unsigned int numLayers; //Num of Layers in the cell after that filter has been applied 
  unsigned int numZtrkSeedCombinations = 0; // Counter for number of seed combinations considered.
  vector<const Stub*> filteredStubs; // Filter Stubs vector to be returned
    
  int FirstSeedLayers[] = {1,2,11,21,3,12,22}; //Allowed layers for the seeding stub

  unsigned int oldNumLay = 0; //Number of Layers counter, used to keep the seed with more layers 

  // Loop over stubs in HT Cell
  for(const Stub* s: stubs){  
    // Create a temporary container for stubs
    vector<const Stub*> tempStubs;
    // Select the first seeding stub
    if(s->psModule() && std::find(std::begin(FirstSeedLayers), std::end(FirstSeedLayers), s->layerId()) != std::end(FirstSeedLayers)){

      numZtrkSeedCombinations++; //Increase cycle counter
      tempStubs.push_back(s); //Push back seed stub in the temporary container
      double sumSeedDist = 0., oldSumSeedDist = 100000.; //Define variable used to estimate the quality of seeds
      // Loop over the remaining stubs in the cell
      for(const Stub* s2: stubs){
	if(s2!=s){
	  // Calculate the correlation factor between the seeding stub s and the considered stub s2
	  double fcorr = 0.;
	  double sum = 0.;
	  for (int i = 0; i < 100; ++i){
	    double zB = -beamWindowZ_ + (0.5+i)*beamWindowZ_/50; 
	    double z1 = s->zTrk() - (chosenRofZFilter_ - s->r())*zB/s->r();
	    double z2 = s2->zTrk() - (chosenRofZFilter_ - s2->r())*zB/s2->r();
	    sum = sum + z1*z2;
	  }
	  sum = (sum/100) - s->zTrk()*s2->zTrk();
	  fcorr = sum/(s->zTrkRes()*s2->zTrkRes());


	  // Check if the zR values of the two stubs (s & s2) are whitin a certain tolerance range defined by strip uncertainty & beam spot length
	  if( fabs(s2->zTrk() - s->zTrk()) < sqrt(s2->zTrkRes()*s2->zTrkRes() + s->zTrkRes()*s->zTrkRes() - fcorr*s->zTrkRes()*s2->zTrkRes() )) {
	    tempStubs.push_back(s2); // Push back s2 if it satisfies the condition
	    sumSeedDist = sumSeedDist + fabs(s2->zTrk() - s->zTrk());  //Increase the seed quality variable
	  }
	}
      }

      sumSeedDist = sumSeedDist/tempStubs.size();

      numLayers = Utility::countLayers(settings_, tempStubs); // Count the number of layers in the temporary stubs container

      // Check if the current seed has more layers then the previous one
      if(numLayers >= oldNumLay){
	// Chech if the current seed has better quality than the previous one
	if(sumSeedDist < oldSumSeedDist){
	  filteredStubs = tempStubs; //Copy the temporary stubs vector in the filteredStubs vector, which will be returned
	  oldSumSeedDist = sumSeedDist; //Update value of oldSumSeedDist
	  oldNumLay = numLayers; //Update value of oldNumLay
	}
      }
    }
  }

  // Print Missing Track information if debug variable is set to 4
  if(settings_->debug()==4){
    std::vector<const Stub* > matchedStubs;
    unsigned int nMatchedLayersBest;
    const TP* matchedTP = Utility::matchingTP(settings_, stubs, nMatchedLayersBest, matchedStubs);
    
    std::vector<const Stub* > matchedFiltStubs;
    unsigned int nMatchedFiltLayersBest;
    const TP* matchedfiltTP = Utility::matchingTP(settings_, filteredStubs, nMatchedFiltLayersBest, matchedFiltStubs);
        
    if(nMatchedLayersBest >= minStubLayers_ && Utility::countLayers(settings_, filteredStubs) < minStubLayers_ ){
      cout << " ******* NOT ENOUGH LAYERS *******" << endl;
      cout << " ====== TP stubs ====== " << endl;
      for(const Stub* st: matchedStubs){
	cout << "z: "<< st->z() << ", r: "<< st->r() << ", id:" << st->layerId() << endl;
      }
      cout << "num layers "<< oldNumLay << endl;
      cout << " ====== Matched TP stubs ====== " << endl;
      for(const Stub* st: filteredStubs){
	cout << "z: "<< st->z() << ", r: "<< st->r() << ", id:" << st->layerId() << endl;
      }
    } else if(nMatchedLayersBest >= minStubLayers_ && nMatchedFiltLayersBest < minStubLayers_ ){
      cout << " ******* NOT ENOUGH MATCHED LAYERS *******" << endl;
      cout << " ====== TP stubs ====== " << endl;
      for(const Stub* st: matchedStubs){
	cout << "z: "<< st->z() << ", r: "<< st->r() << ", id:" << st->layerId() << endl;
      }
      cout << " ====== Matched TP stubs ====== " << endl;
      for(const Stub* st: filteredStubs){
	cout << "z: "<< st->z() << ", r: "<< st->r() << ", id:" << st->layerId() << endl;
      }
    }
  }

  // Note number of seed combinations used for this track.
  numZtrkSeedCombsPerTrk_.push_back(numZtrkSeedCombinations);
   
  return filteredStubs; // Return the filtered stubs vector
}

//=== Produce a filtered collection of stubs on this track candidate that are consistent with a straight line in r-z using tracklet algo.

vector<const Stub*> TrkRZfilter::seedFilter(const std::vector<const Stub*>& stubs, float trkQoverPt) {
  unsigned int numLayers; //Num of Layers in the cell after that filter has been applied 
  std::vector<const Stub*> filtStubs = stubs;   // Copy stubs vector in filtStubs
  bool FirstSeed = true;
  const int FirstSeedLayers[] = {1,2,11,21,3,12,22,4}; //Allowed layers for the first seeding stubs
  const int SecondSeedLayers[] = {1,2,11,3,21,22,12,23,13,4};//Allowed layers for the second seeding stubs
  set<const Stub*> uniqueFilteredStubs;

  unsigned int numSeedCombinations = 0; // Counter for number of seed combinations considered.
  unsigned int numGoodSeedCombinations = 0; // Counter for number of seed combinations considered with z0 within beam spot length.

  vector<const Stub*> filteredStubs; // Filter Stubs vector to be returned
    
  unsigned int oldNumLay = 0; //Number of Layers counter, used to keep the seed with more layers 

  // Loop over stubs in the HT Cell
  for(const Stub* s0: filtStubs){
    // Select the first available seeding stub (r<70)
    if(s0->psModule() && std::find(std::begin(FirstSeedLayers), std::end(FirstSeedLayers), s0->layerId()) != std::end(FirstSeedLayers)) {
      for(const Stub* s1: filtStubs){ 
	if (numGoodSeedCombinations < maxSeedCombinations_) {
	  // Select the second seeding stub (r<90)
	  if(s1->psModule() && s1->layerId() > s0->layerId() &&  std::find(std::begin(SecondSeedLayers), std::end(SecondSeedLayers), s1->layerId()) != std::end(SecondSeedLayers) ){

	    numSeedCombinations++; //Increase filter cycles counter
	
	    double sumSeedDist = 0., oldSumSeedDist = 1000000.; //Define variable used to estimate the quality of seeds
	    vector<const Stub*> tempStubs;  //Create a temporary container for stubs
	    tempStubs.push_back(s0); //Store the first seeding stub in the temporary container
	    tempStubs.push_back(s1); //Store the second seeding stub in the temporary container

	    double z0 = s1->z() + (-s1->z()+s0->z())*s1->r()/(s1->r()-s0->r()); // Estimate a value of z at the beam spot using the two seeding stubs
	    double z0err = s1->zErr() + ( s1->zErr() + s0->zErr() )*s1->r()/fabs(s1->r()-s0->r()) 
	      + fabs(-s1->z()+s0->z())*(s1->rErr()*fabs(s1->r()-s0->r()) + s1->r()*(s1->rErr() + s0->rErr()) )/((s1->r()-s0->r())*(s1->r()-s0->r())); 

            float zTrk = s1->z() + (-s1->z()+s0->z())*(s1->r()-chosenRofZ_)/(s1->r()-s0->r()); // Estimate a value of z at a chosen Radius using the two seeding stubs
            float zTrkErr = s1->zErr() + ( s1->zErr() + s0->zErr() )*fabs(s1->r()-chosenRofZ_)/fabs(s1->r()-s0->r()) 
	      + fabs(-s1->z()+s0->z())*(s1->rErr()*fabs(s1->r()-s0->r()) + fabs(s1->r()-chosenRofZ_)*(s1->rErr() + s0->rErr()) )/((s1->r()-s0->r())*(s1->r()-s0->r()));
        
	    // If z0 is within the beamspot range loop over the other stubs in the cell
	    if (fabs(z0)<=beamWindowZ_+z0err) {
	      // Check track r-z helix parameters are consistent with it being assigned to current rapidity sector (kills duplicates due to overlapping sectors).
	      if ( (! zTrkSectorCheck_) || (zTrk > zTrkMinSector_ - zTrkErr && zTrk < zTrkMaxSector_ + zTrkErr) ) {

		numGoodSeedCombinations++;
		// unsigned int LiD = 0; //Store the layerId of the stub (KEEP JUST ONE STUB PER LAYER)
		// double oldseed = 1000.; //Store the seed value of the current stub (KEEP JUST ONE STUB PER LAYER)

		// Loop over stubs in vector different from the seeding stubs
		for(const Stub* s: filtStubs){
		  if(s!= s0 && s!= s1){
		    // Calculate the seed and its tolerance
		    double seedDist = (s->z() - s1->z())*(s1->r()-s0->r()) - (s->r() - s1->r())*(s1->z() - s0->z());                        
		    double seedDistRes = (s->zErr()+ s1->zErr() )*fabs(s1->r()-s0->r()) + (s->rErr()+s1->rErr())*fabs(s1->z() - s0->z()) 
           		               + (s0->zErr()+s1->zErr())*fabs(s->r() - s1->r()) + (s0->rErr()+s1->rErr())*fabs(s->z() - s1->z());
		    seedDistRes += seedResolution_; // Add extra configurable contribution to assumed resolution.

		    //If seed is lower than the tolerance push back the stub (KEEP JUST ONE STUB PER LAYER, NOT ENABLED BY DEFAULT)
		    // if(fabs(seedDist) <= seedDistRes){
		    //     if(s->layerId()==LiD){
		    //         if(fabs(seedDist)<fabs(oldseed)){
		    //             tempStubs.pop_back();
		    //             tempStubs.push_back(s);
		    //             LiD = s->layerId();
		    //             sumSeedDist = sumSeedDist + fabs(seedDist) - fabs(oldseed);
		    //             oldseed = seedDist;
		    //         }
		    //     } else {
		    //         tempStubs.push_back(s);
		    //         LiD = s->layerId();
		    //         oldseed = seed;
		    //         sumSeedDist = sumSeedDist + fabs(seedDist);
		    //     }               
		    // }
      

		    //If stub lies on the seeding line, store it in the tempstubs vector                          
		    if(fabs(seedDist) <= seedDistRes){
		      tempStubs.push_back(s);
		      sumSeedDist = sumSeedDist + fabs(seedDist); //Increase the seed quality variable
		    }
		  }
		}
	      }
	    }

	    numLayers = Utility::countLayers(settings_, tempStubs); // Count the number of layers in the temporary stubs container
          
	    sumSeedDist = sumSeedDist/(tempStubs.size()); //Measure the average seed quality per stub for the current seed

	    // Check if the current seed has more layers then the previous one (Keep the best seed)
	    if(keepAllSeed_ == false){
	      if(numLayers >= oldNumLay ){
		// Check if the current seed has better quality than the previous one
		if(sumSeedDist < oldSumSeedDist){
		  filteredStubs = tempStubs; //Copy the temporary stubs vector in the filteredStubs vector, which will be returned
		  oldSumSeedDist = sumSeedDist; //Update value of oldSumSeedDist
		  oldNumLay = numLayers; //Update value of oldNumLay
		  estZ0_ = z0; //Store estimated z0
		  estTanLambda_ = (s1->z() -s0->z())/(s1->r()-s0->r()); // Store estimated tanLambda
		  estValid_ = true; 
		}
	      }
	    } else {
	      // Check if the current seed satisfies the minimum layers requirement (Keep all seed algorithm)
	      if (this->trackCandCheck(numLayers, trkQoverPt)) {
		uniqueFilteredStubs.insert(tempStubs.begin(), tempStubs.end()); //Insert the uniqueStub set
		// If these are the first seeding stubs store the values of z0 and tanLambda
		if(FirstSeed){
		  estZ0_ = z0; //Store estimated z0
		  estTanLambda_ = (s1->z() -s0->z())/(s1->r()-s0->r()); // Store estimated tanLambda
		  estValid_ = true;
		  FirstSeed = false; 
		}
	      }
	    }
	  }   
	}
      }
    }
  }
 
  // Copy stubs from the uniqueFilteredStubs set to the filteredStubs vector (Keep all seed algorithm)
  if(keepAllSeed_ == true){
    for (const Stub* stub : uniqueFilteredStubs) {
      filteredStubs.push_back(stub);
    }
  }

  // Print Missing Track information if debug variable is set to 4
  if(settings_->debug()==4){
    std::vector<const Stub* > matchedStubs;
    unsigned int nMatchedLayersBest;
    const TP* matchedTP = Utility::matchingTP(settings_, stubs, nMatchedLayersBest, matchedStubs);
    
    std::vector<const Stub* > matchedFiltStubs;
    unsigned int nMatchedFiltLayersBest;
    const TP* matchedfiltTP = Utility::matchingTP(settings_, filteredStubs, nMatchedFiltLayersBest, matchedFiltStubs);
        
    if(nMatchedLayersBest >= minStubLayers_ && Utility::countLayers(settings_, filteredStubs) < minStubLayers_ ){
      cout << " ******* NOT ENOUGH LAYERS *******" << endl;
      cout << " ====== TP stubs ====== " << endl;
      for(const Stub* st: matchedStubs){
	cout << "z: "<< st->z() << ", r: "<< st->r() << ", id:" << st->layerId() << endl;
      }
      cout << "num layers "<< oldNumLay << endl;
      cout << " ====== Matched TP stubs ====== " << endl;
      for(const Stub* st: filteredStubs){
	cout << "z: "<< st->z() << ", r: "<< st->r() << ", id:" << st->layerId() << endl;
      }
    } else if(nMatchedLayersBest >= minStubLayers_ && nMatchedFiltLayersBest < minStubLayers_ ){
      cout << " ******* NOT ENOUGH MATCHED LAYERS *******" << endl;
      cout << " ====== TP stubs ====== " << endl;
      for(const Stub* st: matchedStubs){
	cout << "z: "<< st->z() << ", r: "<< st->r() << ", id:" << st->layerId() << endl;
      }
      cout << " ====== Matched TP stubs ====== " << endl;
      for(const Stub* st: filteredStubs){
	cout << "z: "<< st->z() << ", r: "<< st->r() << ", id:" << st->layerId() << endl;
      }
    } else if(nMatchedLayersBest >= minStubLayers_ && nMatchedFiltLayersBest >= minStubLayers_ ){
      cout << " ******* Track Found *******" << endl;
      cout << " ====== Cell Stubs ====== " << endl;
      for(const Stub* st: stubs){
	cout << "z: "<< st->digitalStub().iDigi_Z() << ", rT: "<< st->digitalStub().iDigi_Rt() << ", id:" << st->layerId() << endl;
      }
      cout << " ====== Matched TP stubs ====== " << endl;
      for(const Stub* st: filteredStubs){
	cout << "z: "<< st->digitalStub().iDigi_Z() << ", rT: "<< st->digitalStub().iDigi_Rt() << ", id:" << st->layerId() << endl;
      }
    }
  }
    
  // Note number of seed combinations used for this track.
  numSeedCombsPerTrk_.push_back(numSeedCombinations);
  numGoodSeedCombsPerTrk_.push_back(numGoodSeedCombinations);

  return filteredStubs; // Return the filteredStubs vector
}
