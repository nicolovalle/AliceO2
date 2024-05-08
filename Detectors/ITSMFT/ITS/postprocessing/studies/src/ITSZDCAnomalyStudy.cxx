// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "ITSStudies/ITSZDCAnomalyStudy.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DataFormatsParameters/GRPObject.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CCDBTimeStampUtils.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"

#include <set>

#include <TTree.h>

// ZDC
#include "DataFormatsZDC/RecEventFlat.h"


#include "Framework/Task.h"
#include "Framework/Logger.h"

using namespace o2::framework;
using namespace o2::globaltracking;
using GTrackID = o2::dataformats::GlobalTrackID;

namespace o2::its::study
{
class ITSZDCAnomalyStudy : public Task
{
 public:
  ITSZDCAnomalyStudy(std::shared_ptr<DataRequest> dr,
                     std::shared_ptr<o2::base::GRPGeomRequest> gr,
                     bool isMC) : mDataRequest{dr}, mGGCCDBRequest(gr), mUseMC(isMC) {}

  void init(InitContext& ic) final;
  void run(ProcessingContext&) final;
  void endOfStream(EndOfStreamContext&) final;
  void finaliseCCDB(ConcreteDataMatcher&, void*) final;

  // Custom
  void process(o2::globaltracking::RecoContainer& recoData);
  void updateTimeDependentParams(ProcessingContext& pc);

 private:

  void getClusterPatterns(gsl::span<const o2::itsmft::CompClusterExt>&, gsl::span<const unsigned char>&, const o2::itsmft::TopologyDictionary&);
  std::vector<o2::itsmft::ClusterPattern> mPatterns;

  int ChipToLayer(int chip);
  int bcdistance(long orb1, int bc1, long orb2, int bc2);
  std::pair<int,int> findclosestbkg(long orb, int bc, std::map<long,std::set<int>> ZDC);
  
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  std::shared_ptr<DataRequest> mDataRequest;
  bool mUseMC;
  int mStrobe = 594;
  size_t mTFn = 0;
  /*
  std::unique_ptr<TH1F> ZNACall;
  std::unique_ptr<TH1F> ZDCtagBC;
  std::unique_ptr<TH1I> Counters;
  std::unique_ptr<TTree> ITSEvtTree;
  */
  TH1F* ZNACall;
  TH1F* ZDCtagBC;
  TH1I* Counters;
  TTree* ITSEvtTree;
  const o2::itsmft::TopologyDictionary *mDict = nullptr;

  // Tree variables
  int Tbc;
  long Torbit;
  int Tchip;
  int Trofinorbit;
  int TZDCtag;
  int Tclosest_low, Tclosest_up;
  int Tnhit, Tnclus;
  int Tnclus_s20, Tnclus_s100, Tnclus_s150;
  int Tnclus_c20, Tnclus_c100, Tnclus_c128;
  double Tnhit1, Tnhit10;
};

void ITSZDCAnomalyStudy::updateTimeDependentParams(ProcessingContext& pc)
{
  // o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  // static bool initOnceDone = false;
  // if (!initOnceDone) { // this param need to be queried only once
  //   initOnceDone = true;
  //   // mGeom = o2::its::GeometryTGeo::Instance();
  //   // mGeom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot, o2::math_utils::TransformType::T2G));
  // }
}

void ITSZDCAnomalyStudy::init(InitContext& ic)
{
  LOGP(info, "Initializing ITSZDCAnomalyStudy");
  LOGP(info, "Fetching ClusterDictionary");
  auto &mgr = o2::ccdb::BasicCCDBManager::instance();
  mgr.setURL("http://alice-ccdb.cern.ch");
  mgr.setTimestamp(o2::ccdb::getCurrentTimestamp()/*1714900165000*/);
  mDict = mgr.get<o2::itsmft::TopologyDictionary>("ITS/Calib/ClusterDictionary");

  /*
  ZNACall.reset(new TH1F("ZNACall","ZNACall",40,-20,20));
  ZDCtagBC.reset(new TH1F("ZDC tagged BC","ZDC tagged bc",3564,0,3564));
  Counters.reset(new TH1I("Counters","Counters",5,1,6));
  */
  ZNACall = new TH1F("ZNACall","ZNACall",40,-20,20);
  ZDCtagBC = new TH1F("ZDC tagged BC","ZDC tagged bc",3564,0,3564);
  Counters = new TH1I("Counters","Counters",5,1,6);
  
  Counters->GetXaxis()->SetBinLabel(1,"TF"); 
  Counters->GetXaxis()->SetBinLabel(2,"ROF"); 
  Counters->GetXaxis()->SetBinLabel(3,"ZDC evt"); 
  Counters->GetXaxis()->SetBinLabel(4,"ROF-ZDC tagged");
  
  /*
  ITSEvtTree.reset(new TTree("evt","evt"));
  */
  ITSEvtTree = new TTree("evt","evt");

  ITSEvtTree->Branch("orbit",&Torbit,"orbit/L");
  ITSEvtTree->Branch("bc",&Tbc,"bc/I");
  ITSEvtTree->Branch("chip",&Tchip,"chip/I");
  //ITSEvtTree->Branch("rofinorbit",&Trofinorbit,"rofinorbit/I");
  ITSEvtTree->Branch("zdctag",&TZDCtag,"zdctag/I");
  ITSEvtTree->Branch("closest_low",&Tclosest_low,"closest_low/I");
  ITSEvtTree->Branch("closest_up",&Tclosest_up,"closest_up/I");
  ITSEvtTree->Branch("nhit",&Tnhit,"nhit/I");
  ITSEvtTree->Branch("size1",&Tnhit1,"size1/D");
  ITSEvtTree->Branch("size10",&Tnhit10,"size10/D");
  ITSEvtTree->Branch("nclus",&Tnclus,"nclus/I");
  ITSEvtTree->Branch("nclus_s20",&Tnclus_s20,"nclus_s20/I");
  ITSEvtTree->Branch("nclus_s100",&Tnclus_s100,"nclus_s100/I");
  ITSEvtTree->Branch("nclus_s150",&Tnclus_s150,"nclus_s150/I");
  ITSEvtTree->Branch("nclus_c20",&Tnclus_c20,"nclus_c20/I");
  ITSEvtTree->Branch("nclus_c100",&Tnclus_c100,"nclus_c100/I");
  ITSEvtTree->Branch("nclus_c128",&Tnclus_c128,"nclus_c128/I");
  
  
}

  
void ITSZDCAnomalyStudy::endOfStream(EndOfStreamContext&)
{
  LOGP(info, "End of stream for ITSZDCAnomalyStudy");

  std::string outfile1 = "output1.root";
  LOGP(info, "Writing {}", outfile1);
  TFile *F1 = TFile::Open(outfile1.c_str(),"recreate");
  ZNACall->Write();
  ZDCtagBC->Write();
  Counters->Write();
  ITSEvtTree->Write();
  F1->Close();
  delete ZNACall;
  delete ZDCtagBC;
  delete Counters;
  delete ITSEvtTree;

}

void ITSZDCAnomalyStudy::run(ProcessingContext& pc)
{
  LOGP(info, "Running ITSZDCAnomalyStudy on TF: {}", mTFn++);
  o2::globaltracking::RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest.get());
  // updateTimeDependentParams(pc);
  // if (mTFn > 1) return;
  LOGP(info, "Calling process() for TF: {}", mTFn);
  process(recoData);
  
}

void ITSZDCAnomalyStudy::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  return;
}
  /*
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  //   if (matcher == ConcreteDataMatcher("ITS", "CLUSDICT", 0)) {
  //     setClusterDictionary((const o2::itsmft::TopologyDictionary*)obj);
  //     return;
  //   }
}
  */

// Custom area
void ITSZDCAnomalyStudy::process(o2::globaltracking::RecoContainer& recoData)
{
  

  LOGP(info, "Processing RecoContainer");
  Counters->Fill(1);
  LOGP(info, "Retrieving ZDC data");
  auto RecBC = recoData.getZDCBCRecData();
  auto Energy = recoData.getZDCEnergy();
  auto TDCData = recoData.getZDCTDCData();
  auto Info2 = recoData.getZDCInfo();
  LOGP(info, "sizeof ZDC RC: {}, {}, {}, {}", RecBC.size(), Energy.size(), TDCData.size(), Info2.size());

  LOGP(info, "Retrieving ITS clusters");
  auto rofRecVec = recoData.getITSClustersROFRecords();
  auto clusArr = recoData.getITSClusters();
  auto clusPatt = recoData.getITSClustersPatterns();
  LOGP(info, "sizeof ITS RC: {}, {}, {}", clusArr.size(), clusPatt.size(), rofRecVec.size());

  std::map<long,std::set<int>> ZDCtag{}; // ZDCtag[orbit] = <list of bc...>

  // ________________________________________________________________
  // FILLING ZDC ARRAY
  o2::zdc::RecEventFlat ev;
  ev.init(RecBC, Energy, TDCData, Info2);

  int bkgcounter = 0;
  while (ev.next()) {
    
    int32_t itdc = o2::zdc::TDCZNAC; // should be == 0
    
    int nhit = ev.NtdcV(itdc);
    
    for (int32_t ipos = 0; ipos < nhit; ipos++) {
      double mytdc = o2::zdc::FTDCVal * ev.TDCVal[itdc][ipos];
      
      ZNACall->Fill(mytdc);
      
      if (itdc == o2::zdc::TDCZNAC && mytdc > 5.7 && mytdc < 8.7) {

	// Backgroud event found here!
	bkgcounter++;
	Counters->Fill(3);
	ZDCtagBC->Fill(ev.ir.bc);
	long zdcorbit = (long)ev.ir.orbit;

	if (ZDCtag.find(zdcorbit) != ZDCtag.end()){
	  bool double_count_bkg = ZDCtag[zdcorbit].insert((int)ev.ir.bc).second;	  
	  if (double_count_bkg){
	    LOGP(warning,"Multiple ZDC counts in the same orbit/bc {}/{}",zdcorbit,ev.ir.bc);
	  }
	}
	else {
	  std::set<int> zdcbcs{(int)ev.ir.bc};
	  ZDCtag[zdcorbit] = zdcbcs;
	}
	
      }
    }
  } // end of while ev.next()

  LOGP(info,"Found {} background envents from ZNAC",bkgcounter++);
  //__________________________________________________________________



  std::map<long,bool> ChipSeen{}; // ChipSeen[orbit*3546 + bc] = true/false
  /*
  for (auto& rofRec : rofRecVec){
    int l1bc = (int)rofRec.getBCData().bc;
    long l1orbit = (long)rofRec.getBCData().orbit;
    auto cl
  */

  
  getClusterPatterns(clusArr, clusPatt, *mDict);
  //auto pattIt = clusPatt.begin();

  bool printFirstOrbit = false;
  
  for (auto& rofRec : rofRecVec){  // size: ROFs in TF

    Counters->Fill(2);

    auto clustersInRof = rofRec.getROFData(clusArr);
    auto patternsInRof = rofRec.getROFData(mPatterns); 

    Tbc = (int)rofRec.getBCData().bc;
    Torbit = (long)rofRec.getBCData().orbit;

    if (!printFirstOrbit){
      LOGP(info,"First of TF: ITS orbit/bc {}/{}",Torbit,Tbc);
      printFirstOrbit = true;
    }

    Trofinorbit = (int)Tbc/mStrobe;

    TZDCtag = 0;
    if (ZDCtag.find((long)Torbit) != ZDCtag.end()){
      for (auto zbc : ZDCtag[(long)Torbit]){
	if ((int)(zbc/mStrobe) == Trofinorbit){
	  TZDCtag = 1;
	  Counters->Fill(4);
	  break;
	}
      }
    }

    std::pair<int,int> closestbkg = findclosestbkg(Torbit, Tbc, ZDCtag);

    Tclosest_low = closestbkg.first;
    Tclosest_up = closestbkg.second;

    if (TZDCtag > 0) LOGP(info,"Closest bkg: asking {}/{}, returning {} , {}. ZDCTag {}",Torbit,Tbc,Tclosest_low,Tclosest_up,TZDCtag);
    //else LOGP(info,"DEBUGcl bkg: asking {}/{}, returning {} , {}. ZDCTag {}",Torbit,Tbc,Tclosest_low,Tclosest_up,TZDCtag);

    std::set<int> AvailableChips{};
    std::map<int, std::vector<int>> MAPsize{};
    std::map<int, std::vector<int>> MAPcols{};

    for (int iclus = 0; iclus < clustersInRof.size(); iclus++){

      const auto& compClus = clustersInRof[iclus];

      auto chipid = compClus.getSensorID();

      if (ChipToLayer(chipid) > 2){
	continue;
      }

      int npix = 0;
      int colspan = 0;

      // Cluster dictionary ////
      /*
      auto pattID = compClus.getPatternID();
      if (pattID == o2::itsmft::CompCluster::InvalidPatternID || mDict->isGroup(pattID)){
	o2::itsmft::ClusterPattern patt(pattIt);
	npix = patt.getNPixels();
	colspan = patt.getColumnSpan();
      }
      else{
	npix = mDict->getNpixels(pattID);
	o2::itsmft::ClusterPattern patter = mDict->getPattern(pattID);
	colspan = patter.getColumnSpan();
      }
      */
      //////////////////
      // Otherwise:
      auto patti = patternsInRof[iclus];
      npix = patti.getNPixels();
      colspan = patti.getColumnSpan();
      /////////////////////////////////////////


      
	
      
     

      bool newchip = AvailableChips.insert(chipid).second;
      if (newchip) {
	MAPsize[chipid] = std::vector<int>{};
	MAPcols[chipid] = std::vector<int>{};
      }

      MAPsize[chipid].push_back(npix);
      MAPcols[chipid].push_back(colspan);
      
 
    } // end of loop over clusters in rof

    for (int ic : AvailableChips){

      Tchip = ic;

      Tnclus = MAPsize[ic].size();

      std::sort(MAPsize[ic].begin(), MAPsize[ic].end(), std::greater<int>());
      //std::sort(MAPcols[ic].begin(), MAPcols[ic].end(), std::greater<int>());

      Tnhit = Tnclus_s20 = Tnclus_s100 = Tnclus_s150 = 0;
      Tnhit1 = Tnhit10 = 0.;
      Tnclus_c20 = Tnclus_c100 = Tnclus_c128 = 0;
      
      int nclus10 = 0, nclus1 = 0;
      
      for (int nh: MAPsize[ic]){

	Tnhit += nh;
	
	if (nclus10 < 10){
	  nclus10++;
	  Tnhit10 += 1.*nh;
	}

	if (nclus1 < 1){
	  nclus1++;
	  Tnhit1 += 1.*nh;
	}

	Tnclus_s20  += (nh >= 20);
	Tnclus_s100 += (nh >= 100);
	Tnclus_s150 += (nh >= 150);	
      }

      Tnhit10 = (nclus10 == 0) ? 0. : 1.*Tnhit10/nclus10;

      for (int nc: MAPcols[ic]){
	Tnclus_c20  += (nc >= 20);
	Tnclus_c100 += (nc >= 100);
	Tnclus_c128 += (nc >= 128);	
      }


      // filling the Event Tree
      ITSEvtTree->Fill();

    } // end of loop over available chips


  
  } // end of loop over ROFs
    
  
    
  

  // plot interessanti:
  // Cluster size (averged or max) vs chip
  // Cluster length vs chip
  // Cluster prima d un evento in cui il chip è vuoto
  // numero di chip con cluster grande
  //    (x2 --> tutti gli eventi o con solo gli eventi ZDC-tagged)

  // Osservabile vs DeltaT (evento ITS e evento ZDC-tagged più vicino).

  // print tagged orbit-bc
  if (ZDCtag.size()>0){
    for (auto ob : ZDCtag){
      for (auto b : ob.second){
	LOGP(info,"ZNAC background on orbit/bc {}/{}",ob.first,b);
      }
    }
  }
}

int ITSZDCAnomalyStudy::ChipToLayer(int chip){
   if (chip < 108)
    return 0;
  if (chip < 252)
    return 1;
  if (chip < 432)
    return 2;
  if (chip < 3120)
    return 3;
  if (chip < 6480)
    return 4;
  if (chip < 14712)
    return 5;
  return 6;
}
  
int ITSZDCAnomalyStudy::bcdistance(long orb1, int bc1, long orb2, int bc2){ // returns "2" - "1"
  if (orb1 == orb2) return bc2-bc1;
  long dOrb = orb2-orb1;
  int dOBC = static_cast<int>(dOrb * 3564);
  if (dOrb > 0) return dOBC - bc1 + bc2;
  else return 0-ITSZDCAnomalyStudy::bcdistance(orb2,bc2,orb1,bc1);
}


std::pair<int, int> ITSZDCAnomalyStudy::findclosestbkg(long orb, int bc, std::map<long,std::set<int>> ZDC)
  {
    if (ZDC.size() == 0){
      return std::make_pair(999999,-999999);
    }
    std::vector<int> distances{-999999};
    //std::vector<long> debug_orbit{-999999};
    //std::vector<int> debug_bc{-999999};
    for (auto ob : ZDC){
      for (auto b : ob.second){
	distances.push_back(bcdistance(orb,bc,ob.first,b));
	//debug_orbit.push_back(ob.first);
	//debug_bc.push_back(b);
      }
    }

    distances.push_back(999999);
    //debug_orbit.push_back(-999999);
    //debug_bc.push_back(-999999);

    int toret1, toret2;

    int i=0;
    for (i=0; i<distances.size()-1; i++){
      if (distances[i] < 0 && distances[i+1]>=0){
	toret1 = distances[i];
	toret2 = distances[i+1];	
	break;
      }
    }

    //LOGP(info,"Closest bkg function. Asked {}/{} , returning distances {}({}/{}) , {}({}/{})", orb, bc, toret1, debug_orbit[i], debug_bc[i], toret2, debug_orbit[i+1], debug_bc[i+1]);

    if (toret1 == -999999) toret1 = 999999;
    if (toret2 == 999999) toret2 = -999999;

    return std::make_pair(toret1,toret2);
  }


void ITSZDCAnomalyStudy::getClusterPatterns(gsl::span<const o2::itsmft::CompClusterExt>& ITSclus, gsl::span<const unsigned char>& ITSpatt, const o2::itsmft::TopologyDictionary& mdict)
{
  mPatterns.clear();
  mPatterns.reserve(ITSclus.size());
  auto pattIt = ITSpatt.begin();

  for (unsigned int iClus{0}; iClus < ITSclus.size(); ++iClus) {
    auto& clus = ITSclus[iClus];

    auto pattID = clus.getPatternID();
    o2::itsmft::ClusterPattern patt;

    if (pattID == o2::itsmft::CompCluster::InvalidPatternID || mdict.isGroup(pattID)) {
      patt.acquirePattern(pattIt);
    } else {
      patt = mdict.getPattern(pattID);
    }

    mPatterns.push_back(patt);
  }
}





// getter
  DataProcessorSpec getITSZDCAnomalyStudy(mask_t srcTracksMask, mask_t srcClustersMask, bool useMC)
{

  // std::cout<<"DEBBUG track and clus masks "<<srcTracksMask<<" "<<srcClustersMask<<" is ZDC in tracks: "<<(srcTracksMask & GTrackID::getSourcesMask("ZDC"))<<" is ITS in clus: "<<(srcClustersMask & GTrackID::getSourcesMask("ITS"))<<std::endl;
  
  std::vector<OutputSpec> outputs;
  auto dataRequest = std::make_shared<DataRequest>();
  dataRequest->requestClusters(srcClustersMask, useMC); // working
  //dataRequest->requestTracks(GTrackID::getSourcesMask("ZDC"), useMC); //working
  
 
  dataRequest->requestTracks(srcTracksMask, useMC);

  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                              true,                              // GRPECS=true
                                                              false,                             // GRPLHCIF
                                                              false,                             // GRPMagField
                                                              false,                             // askMatLUT
                                                              o2::base::GRPGeomRequest::Aligned, // geometry
                                                              dataRequest->inputs,
                                                              true);
  return DataProcessorSpec{
    "its-zdc-anomaly-study",
    dataRequest->inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<ITSZDCAnomalyStudy>(dataRequest, ggRequest, useMC)},
    Options{}};
}

} // namespace o2::its::study
