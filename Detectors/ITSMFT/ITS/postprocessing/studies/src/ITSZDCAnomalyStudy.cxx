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
#include <algorithm>
#include <limits>

#include <TTree.h>
#include <TH2.h>

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
  void save_and_reset();

  // Custom
  void process(o2::globaltracking::RecoContainer& recoData);
  void updateTimeDependentParams(ProcessingContext& pc);

 private:
  void getClusterPatterns(gsl::span<const o2::itsmft::CompClusterExt>&, gsl::span<const unsigned char>&, const o2::itsmft::TopologyDictionary&);
  std::vector<o2::itsmft::ClusterPattern> mPatterns;

  int ChipToLayer(int chip);
  double ChipToPhi(int chip);

  int NStaves[7] = {12, 16, 20, 24, 30, 42, 48};
  int N_CHIP_IB = 9 * (NStaves[0] + NStaves[1] + NStaves[2]);
  int bcdistance(long orb1, int bc1, long orb2, int bc2);
  std::pair<long, int> shiftbc(long orb, int bc, int shift);
  std::pair<int, int> findclosestbkg(long orb, int bc, std::map<long, std::set<int>> ZDC);

  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  std::shared_ptr<DataRequest> mDataRequest;
  bool mUseMC;
  int mStrobeFallBack = 594;
  int mStrobe = mStrobeFallBack;
  int mTFn = 0;
  int mTF_first_after_dump = 1;

  // custom settings // to make configurables...
  std::pair<double, double> TimeWindowZDC = std::make_pair(5., 11.); // used up to 27 oct 24
  std::pair<double, double> TimeWindowZNAr = std::make_pair(5.64507, 9.64507);
  std::pair<double, double> TimeWindowZNCr = std::make_pair(4.21299, 8.21299);
  std::pair<double, double> TimeWindowZNAl = std::make_pair(-11.3549, -8.35493);
  std::pair<double, double> TimeWindowZNCl = std::make_pair(-12.787, -9.78701);

  int mDumpEveryTF = 1000; // use -1 to save only at the end
  int mSkimmedOnlyAfterTF = 15000;
  std::string mOutput1  = "chipevents";
  std::string mOutput11 = "chipeventstarget";
  std::string mOutput2  = "staveevents";
  std::string mOutput3  = "rofevents";

  /*
  std::unique_ptr<TH1F> ZNACall;
  std::unique_ptr<TH1F> ZDCtagBC;
  std::unique_ptr<TH1I> Counters;
  std::unique_ptr<TTree> ITSChipEvtTree;
  */
  TH1F* TimeWindowCut;
  TH1F* ZNACall;
  TH1F* ZNCCall;
  TH1F* ZDCAtagBC;
  TH1F* ZDCCtagBC;
  TH1I* Counters;
  TH2D* ClusterShape;
  TH2D* ClusterShapeTO;
  TTree* ITSChipEvtTree;
  TTree* ITSChipEvtTargetTree;
  TTree* ITSStaveEvtTree;
  TTree* ITSROFEvtTree;
  const o2::itsmft::TopologyDictionary* mDict = nullptr;

  // Tree variables
  int Tbc;
  long Torbit;
  int Tchip;
  double Tphi;
  int TZDCtag;
  int TZDCtag_extended; // looks also into previous and next ROF
  int Tclosest_low, Tclosest_up;
  int Tnhit, Tnclus, Tnhit_no1pix, Tnclus_no1pix;
  double Tstdhit, Tstdhit_no1pix;
  int Tnclus_s20, Tnclus_s100, Tnclus_s150;
  int Tnclus_c20, Tnclus_c100, Tnclus_c128;
  int Tnclus_target;
  double Tnhit1, Tnhit10;
  int Tmissingafter, Tmissingafter2;
  int Tmincol, Tmaxcol;

  // Stave tree variables
  int Sstave;
  int Snchip;
  double Sphi;
  int Snhit, Snclus;
  int Snclus_s20, Snclus_s100, Snclus_s150;
  int Snclus_c20, Snclus_c100, Snclus_c128, Snclus_target;
  double Snhit1, Snhit10;

  // ROF tree variables
  int Ritstag, RitstagTO;
  int Rnchip, RnchipTO;
  int RnchipafterTO, Rnchipafter;
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
  auto& mgr = o2::ccdb::BasicCCDBManager::instance();
  mgr.setURL("http://alice-ccdb.cern.ch");
  mgr.setTimestamp(o2::ccdb::getCurrentTimestamp() /*1714900165000*/);
  mDict = mgr.get<o2::itsmft::TopologyDictionary>("ITS/Calib/ClusterDictionary");

  LOGP(info, "Setting up trees and histograms. They will be dumped on file every {} TFs", mDumpEveryTF);
  /*
  ZNACall.reset(new TH1F("ZNACall","ZNACall",40,-20,20));
  ZDCtagBC.reset(new TH1F("ZDC tagged BC","ZDC tagged bc",3564,0,3564));
  Counters.reset(new TH1I("Counters","Counters",5,1,6));
  */
  TimeWindowCut = new TH1F("ZDC bkg region", "ZNAr, ZNCr, -ZNAl, -ZNCl", 8, 0, 8);
  TimeWindowCut->SetBinContent(1, TimeWindowZNAr.first);
  TimeWindowCut->SetBinContent(2, TimeWindowZNAr.second);
  TimeWindowCut->SetBinContent(3, TimeWindowZNCr.first);
  TimeWindowCut->SetBinContent(4, TimeWindowZNCr.second);
  TimeWindowCut->SetBinContent(5, -TimeWindowZNAl.first);
  TimeWindowCut->SetBinContent(6, -TimeWindowZNAl.second);
  TimeWindowCut->SetBinContent(7, -TimeWindowZNCl.first);
  TimeWindowCut->SetBinContent(8, -TimeWindowZNCl.second);
  ZNACall = new TH1F("ZNACall", "ZNACall", 40, -20, 20);
  ZNCCall = new TH1F("ZNCCall", "ZNCCall", 40, -20, 20);
  ZDCAtagBC = new TH1F("ZDCA tagged BC", "ZDCA tagged bc", 3564, 0, 3564);
  ZDCCtagBC = new TH1F("ZDCC tagged BC", "ZDCC tagged bc", 3564, 0, 3564);
  Counters = new TH1I("Counters", "Counters", 20, 1, 21);
  ClusterShape = new TH2D("ClusterShape", "Col and Row span of cluster npix > 1; colspan; rowspan", 130, 0, 130, 130, 0, 130);
  ClusterShapeTO = new TH2D("ClusterShapeTO", "Col and Row span of cluster in IB-TO, npix > 1; colspan; rowspan", 130, 0, 130, 130, 0, 130);

  Counters->GetXaxis()->SetBinLabel(1, "TF");
  Counters->GetXaxis()->SetBinLabel(2, "ROF");
  Counters->GetXaxis()->SetBinLabel(3, "ZDCA evt");
  Counters->GetXaxis()->SetBinLabel(4, "ROF-ZDCA tagged");
  Counters->GetXaxis()->SetBinLabel(5, "ITStag any");
  Counters->GetXaxis()->SetBinLabel(6, "ITStag TO");
  Counters->GetXaxis()->SetBinLabel(7, "ITStag any + ZDC");
  Counters->GetXaxis()->SetBinLabel(8, "ITStag TO + ZDC");
  Counters->GetXaxis()->SetBinLabel(9, "ZDCC evt");
  Counters->GetXaxis()->SetBinLabel(10, "ROF-ZDCC tagged");


  ITSChipEvtTree = new TTree("chipevt", "chipevt");
  ITSStaveEvtTree = new TTree("staveevt", "staveevt");
  ITSROFEvtTree = new TTree("rofevt", "rofevt");

  // Chip event branches
  ITSChipEvtTree->Branch("TFprogress",&mTFn,"TFprogress/I");
  ITSChipEvtTree->Branch("orbit", &Torbit, "orbit/L");
  ITSChipEvtTree->Branch("bc", &Tbc, "bc/I");
  ITSChipEvtTree->Branch("chip", &Tchip, "chip/I");
  ITSChipEvtTree->Branch("phi", &Tphi, "phi/D");
  ITSChipEvtTree->Branch("zdctag", &TZDCtag, "zdctag/I");
  // ITSChipEvtTree->Branch("zdctag_ext",&TZDCtag_extended,"zdctag/I");
  // ITSChipEvtTree->Branch("closest_low",&Tclosest_low,"closest_low/I");
  // ITSChipEvtTree->Branch("closest_up",&Tclosest_up,"closest_up/I");
  ITSChipEvtTree->Branch("nhit", &Tnhit, "nhit/I");
  ITSChipEvtTree->Branch("nhit_no1pix", &Tnhit_no1pix, "nhit_no1pix/I");
  // ITSChipEvtTree->Branch("stdhit",&Tstdhit,"stdhit/D");
  // ITSChipEvtTree->Branch("stdhit_no1pix",&Tstdhit_no1pix,"stdhit_no1pix/D");
  ITSChipEvtTree->Branch("size1", &Tnhit1, "size1/D");
  ITSChipEvtTree->Branch("size10", &Tnhit10, "size10/D");
  ITSChipEvtTree->Branch("nclus", &Tnclus, "nclus/I");
  ITSChipEvtTree->Branch("nclus_no1pix", &Tnclus_no1pix, "nclus_no1pix/I");
  // ITSChipEvtTree->Branch("nclus_s20",&Tnclus_s20,"nclus_s20/I");
  // ITSChipEvtTree->Branch("nclus_s100",&Tnclus_s100,"nclus_s100/I");
  // ITSChipEvtTree->Branch("nclus_s150",&Tnclus_s150,"nclus_s150/I");
  // ITSChipEvtTree->Branch("nclus_c20",&Tnclus_c20,"nclus_c20/I");
  // ITSChipEvtTree->Branch("nclus_c100",&Tnclus_c100,"nclus_c100/I");
  // ITSChipEvtTree->Branch("nclus_c128",&Tnclus_c128,"nclus_c128/I");
  ITSChipEvtTree->Branch("nclus_target", &Tnclus_target, "nclus_target/I");
  ITSChipEvtTree->Branch("missingafter", &Tmissingafter, "missingafter/I");
  ITSChipEvtTree->Branch("missingafter2", &Tmissingafter2, "missingafter2/I");
  ITSChipEvtTree->Branch("mincol", &Tmincol, "mincol/I");
  ITSChipEvtTree->Branch("maxcol", &Tmaxcol, "maxcol/I");
  
  ITSChipEvtTargetTree = ITSChipEvtTree->CloneTree(0);
  ITSChipEvtTargetTree->SetName("chipevttarget");
  ITSChipEvtTargetTree->SetTitle("chipevttarget");


  ITSStaveEvtTree->Branch("orbit", &Torbit, "orbit/L");
  ITSStaveEvtTree->Branch("bc", &Tbc, "bc/I");
  ITSStaveEvtTree->Branch("stave", &Sstave, "stave/I");
  ITSStaveEvtTree->Branch("phi", &Sphi, "phi/D");
  ITSStaveEvtTree->Branch("zdctag", &TZDCtag, "zdctag/I");
  ITSStaveEvtTree->Branch("zdctag_ext", &TZDCtag_extended, "zdctag/I");
  // ITSStaveEvtTree->Branch("closest_low",&Tclosest_low,"closest_low/I");
  // ITSStaveEvtTree->Branch("closest_up",&Tclosest_up,"closest_up/I");
  ITSStaveEvtTree->Branch("nchip", &Snchip, "nchip/I");
  ITSStaveEvtTree->Branch("nhit", &Snhit, "nhit/I");
  ITSStaveEvtTree->Branch("size1", &Snhit1, "size1/D");
  ITSStaveEvtTree->Branch("size10", &Snhit10, "size10/D");
  ITSStaveEvtTree->Branch("nclus", &Snclus, "nclus/I");
  // ITSStaveEvtTree->Branch("nclus_s20",&Snclus_s20,"nclus_s20/I");
  // ITSStaveEvtTree->Branch("nclus_s100",&Snclus_s100,"nclus_s100/I");
  // ITSStaveEvtTree->Branch("nclus_s150",&Snclus_s150,"nclus_s150/I");
  // ITSStaveEvtTree->Branch("nclus_c20",&Snclus_c20,"nclus_c20/I");
  // ITSStaveEvtTree->Branch("nclus_c100",&Snclus_c100,"nclus_c100/I");
  ITSStaveEvtTree->Branch("nclus_target", &Snclus_target, "nclus_target/I");

  ITSROFEvtTree->Branch("orbit", &Torbit, "orbit/L");
  ITSROFEvtTree->Branch("bc", &Tbc, "bc/I");
  ITSROFEvtTree->Branch("zdctag", &TZDCtag, "zdctag/I");
  ITSROFEvtTree->Branch("zdctag_ext", &TZDCtag_extended, "zdctag/I");
  ITSROFEvtTree->Branch("closest_low", &Tclosest_low, "closest_low/I");
  ITSROFEvtTree->Branch("closest_up", &Tclosest_up, "closest_up/I");
  ITSROFEvtTree->Branch("itstag", &Ritstag, "itstag/I");
  //ITSROFEvtTree->Branch("itstagTO", &RitstagTO, "itstagTO/I");
  ITSROFEvtTree->Branch("nchip", &Rnchip, "nchip/I");
  //ITSROFEvtTree->Branch("nchipTO", &RnchipTO, "nchipTO/I");
  ITSROFEvtTree->Branch("nchipafter", &Rnchipafter, "nchipafter/I");
}

void ITSZDCAnomalyStudy::save_and_reset()
{

  std::string outfile11 = mOutput11 + "_" + std::to_string(mTF_first_after_dump) + "_" + std::to_string(mTFn) + ".root";
  LOGP(info, "Writing ROOT file {}", outfile11);
  TFile* F11 = TFile::Open(outfile11.c_str(), "recreate");
  TimeWindowCut->Write();
  ZNACall->Write();
  ZNCCall->Write();
  ZDCAtagBC->Write();
  ZDCCtagBC->Write();
  Counters->Write();
  ClusterShape->Write();
  ClusterShapeTO->Write();
  ITSChipEvtTargetTree->Write();
  F11->Close();
  delete F11;

  if (mTFn <= mSkimmedOnlyAfterTF){
    std::string outfile1 = mOutput1 + "_" + std::to_string(mTF_first_after_dump) + "_" + std::to_string(mTFn) + ".root";
    LOGP(info, "Writing ROOT file {}", outfile1);
    TFile* F1 = TFile::Open(outfile1.c_str(), "recreate");
    TimeWindowCut->Write();
    ZNACall->Write();
    ZNCCall->Write();
    ZDCAtagBC->Write();
    ZDCCtagBC->Write();
    Counters->Write();
    ClusterShape->Write();
    ClusterShapeTO->Write();
    ITSChipEvtTree->Write();  // chip events and the skimmed one
    ITSChipEvtTargetTree->Write();
    F1->Close();
    delete F1;
    
    std::string outfile2 = mOutput2 + "_" + std::to_string(mTF_first_after_dump) + "_" + std::to_string(mTFn) + ".root";
    LOGP(info, "Writing ROOT file {}", outfile2);
    TFile* F2 = TFile::Open(outfile2.c_str(), "recreate");
    TimeWindowCut->Write();
    ZNACall->Write();
    ZNCCall->Write();
    ZDCAtagBC->Write();
    ZDCCtagBC->Write();
    Counters->Write();
    ClusterShape->Write();
    ClusterShapeTO->Write();
    ITSStaveEvtTree->Write();  // stave event
    F2->Close();
    delete F2;

    std::string outfile3 = mOutput3 + "_" + std::to_string(mTF_first_after_dump) + "_" + std::to_string(mTFn) + ".root";
    LOGP(info, "Writing ROOT file {}", outfile3);
    TFile* F3 = TFile::Open(outfile3.c_str(), "recreate");
    TimeWindowCut->Write();
    ZNACall->Write();
    ZNCCall->Write();
    ZDCAtagBC->Write();
    ZDCCtagBC->Write();
    Counters->Write();
    ClusterShape->Write();
    ClusterShapeTO->Write();
    ITSROFEvtTree->Write();  // rof event
    F3->Close();
    delete F3;
  }

  LOGP(info, "Resetting historgrams and trees");
  // Delete clears data but keep the branch setup intact
  ITSChipEvtTree->Reset(); // Delete("");
  ITSChipEvtTargetTree->Reset(); // Delete("");
  ITSROFEvtTree->Reset(); // Delete("");
  ITSStaveEvtTree->Reset(); // Delete("");
  TimeWindowCut->Reset();
  ZNACall->Reset();
  ZNCCall->Reset();
  ZDCAtagBC->Reset();
  ZDCCtagBC->Reset();
  Counters->Reset();
  ClusterShape->Reset();
  ClusterShapeTO->Reset();

  mTF_first_after_dump = mTFn + 1;
}

void ITSZDCAnomalyStudy::endOfStream(EndOfStreamContext&)
{
  LOGP(info, "End of stream for ITSZDCAnomalyStudy");

  save_and_reset();

  delete TimeWindowCut;
  delete ZNACall;
  delete ZNCCall;
  delete ZDCAtagBC;
  delete ZDCCtagBC;
  delete Counters;
  delete ITSChipEvtTree;
  delete ITSChipEvtTargetTree;
  delete ITSStaveEvtTree;
  delete ITSROFEvtTree;
}

void ITSZDCAnomalyStudy::run(ProcessingContext& pc)
{

  if (mTFn == std::numeric_limits<int>::max()){
    LOGP(error, "Refusing to prcess more than {} TFs. Skipping all next events", mTFn);
    return;
  }
  
  mTFn++;

  if (mDumpEveryTF > 0 && mTFn > 0 && (mTFn % mDumpEveryTF) == 0) {
    LOGP(info, "Reached TF #{}. Exporting new root files", mTFn);
    save_and_reset();
  }

  //LOGP(info, "Running ITSZDCAnomalyStudy on TF: {}", mTFn);
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

  if (rofRecVec.size() == 576 || rofRecVec.size() == 192) {
    mStrobe = 3564 / (rofRecVec.size() / 32);
    LOGP(info, "Assuimg TF length = 32 orbits and setting strobe length to {} bc", mStrobe);
  } else {
    mStrobe = mStrobeFallBack;
    LOGP(warning, "Unforeseen number of ROFs in the loop. Using the strobe length fall back value {}", mStrobe);
  }

  std::map<long, std::set<int>> ZNArtag{}; // ZDCAtag[orbit] = <list of bc...>
  std::map<long, std::set<int>> ZNCrtag{};
  std::map<long, std::set<int>> ZNAltag{};
  std::map<long, std::set<int>> ZNCltag{};

  // ________________________________________________________________
  // FILLING ZDC ARRAY
  o2::zdc::RecEventFlat ev;
  ev.init(RecBC, Energy, TDCData, Info2);

  int bkgcounterAr = 0, bkgcounterCr = 0;
  int bkgcounterAl = 0, bkgcounterCl = 0;
  while (ev.next()) {

    int32_t itdcA = o2::zdc::TDCZNAC; // should be == 0
    int32_t itdcC = o2::zdc::TDCZNCC;
    long zdcorbit = (long)ev.ir.orbit;

    // ZDC - A side
    int nhitA = ev.NtdcV(itdcA);
    for (int32_t ipos = 0; ipos < nhitA; ipos++) {

      double mytdc = o2::zdc::FTDCVal * ev.TDCVal[itdcA][ipos];

      ZNACall->Fill(mytdc);

      if (mytdc >= TimeWindowZNAr.first && mytdc <= TimeWindowZNAr.second) { // original from Pietro: [5.7 - 8.7]

        // Backgroud event found here!
        bkgcounterAr++;
        Counters->Fill(3);
        ZDCAtagBC->Fill(ev.ir.bc);

        if (ZNArtag.find(zdcorbit) != ZNArtag.end()) {
          bool double_count_bkg = ZNArtag[zdcorbit].insert((int)ev.ir.bc).second;
          if (double_count_bkg) {
            LOGP(warning, "Multiple ZDCAr counts in the same orbit/bc {}/{}", zdcorbit, ev.ir.bc);
          }
        } else {
          std::set<int> zdcbcs{(int)ev.ir.bc};
          ZNArtag[zdcorbit] = zdcbcs;
        }

      } // and of ZNAr time window

      if (mytdc >= TimeWindowZNAl.first && mytdc <= TimeWindowZNAl.second) {

        // Backgroud event found here!
        bkgcounterAl++;
        Counters->Fill(3);
        ZDCAtagBC->Fill(ev.ir.bc);

        if (ZNAltag.find(zdcorbit) != ZNAltag.end()) {
          bool double_count_bkg = ZNAltag[zdcorbit].insert((int)ev.ir.bc).second;
          if (double_count_bkg) {
            LOGP(warning, "Multiple ZDCAl counts in the same orbit/bc {}/{}", zdcorbit, ev.ir.bc);
          }
        } else {
          std::set<int> zdcbcs{(int)ev.ir.bc};
          ZNAltag[zdcorbit] = zdcbcs;
        }

      } // and of ZNAl time window
    }

    // ZDC - C side
    int nhitC = ev.NtdcV(itdcC);
    for (int32_t ipos = 0; ipos < nhitC; ipos++) {

      double mytdc = o2::zdc::FTDCVal * ev.TDCVal[itdcC][ipos];

      ZNCCall->Fill(mytdc);

      if (mytdc >= TimeWindowZNCr.first && mytdc <= TimeWindowZNCr.second) {

        // Backgroud event found here!
        bkgcounterCr++;
        Counters->Fill(9);
        ZDCCtagBC->Fill(ev.ir.bc);

        if (ZNCrtag.find(zdcorbit) != ZNCrtag.end()) {
          bool double_count_bkg = ZNCrtag[zdcorbit].insert((int)ev.ir.bc).second;
          if (double_count_bkg) {
            LOGP(warning, "Multiple ZNCr counts in the same orbit/bc {}/{}", zdcorbit, ev.ir.bc);
          }
        } else {
          std::set<int> zdcbcs{(int)ev.ir.bc};
          ZNCrtag[zdcorbit] = zdcbcs;
        }

      } // end of ZNCr time window

      if (mytdc >= TimeWindowZNCl.first && mytdc <= TimeWindowZNCl.second) {

        // Backgroud event found here!
        bkgcounterCl++;
        Counters->Fill(9);
        ZDCCtagBC->Fill(ev.ir.bc);

        if (ZNCltag.find(zdcorbit) != ZNCltag.end()) {
          bool double_count_bkg = ZNCltag[zdcorbit].insert((int)ev.ir.bc).second;
          if (double_count_bkg) {
            LOGP(warning, "Multiple ZNCl counts in the same orbit/bc {}/{}", zdcorbit, ev.ir.bc);
          }
        } else {
          std::set<int> zdcbcs{(int)ev.ir.bc};
          ZNCltag[zdcorbit] = zdcbcs;
        }

      } // end of ZNCl time window
    }

  } // end of while ev.next()

  LOGP(info, "Found background envents from ZNAright/left {}/{} -- from ZNCright/left {}/{}", bkgcounterAr, bkgcounterAl, bkgcounterCr, bkgcounterCl);
  //__________________________________________________________________

  getClusterPatterns(clusArr, clusPatt, *mDict);
  // auto pattIt = clusPatt.begin();

  int inTFROFcounter = -1;

  std::vector<bool> ChipSeenInThisROF(N_CHIP_IB, false);  // ChipSeenInThisROF[chipid] = true/false
  std::vector<bool> ChipSeenInLastROF(N_CHIP_IB, false);  // ChipSeenInLastROF[chipid] = true/false
  std::vector<bool> ChipSeenInLast2ROF(N_CHIP_IB, false); // ChipSeenInLast2ROF[chipid] = true/false

  // for (auto& rofRec : rofRecVec){  // size: ROFs in TF
  for (auto it = rofRecVec.rbegin(); it != rofRecVec.rend(); ++it) {

    auto& rofRec = *it;

    inTFROFcounter++;

    Counters->Fill(2);

    ChipSeenInLast2ROF = ChipSeenInLastROF;
    ChipSeenInLastROF = ChipSeenInThisROF;
    std::fill(ChipSeenInThisROF.begin(), ChipSeenInThisROF.end(), false);

    auto clustersInRof = rofRec.getROFData(clusArr);
    auto patternsInRof = rofRec.getROFData(mPatterns);

    Tbc = (int)rofRec.getBCData().bc;
    Torbit = (long)rofRec.getBCData().orbit;

    if (inTFROFcounter < 1) {
      LOGP(info, "First of TF: ITS orbit/bc {}/{}", Torbit, Tbc);
    }

    int eff_bc = Tbc + 60;
    long eff_orbit = Torbit;

    if (eff_bc > 3564) {
      eff_bc -= 3564;
      eff_orbit += 1;
    }

    bool isZNArtagged = false;
    if (ZNArtag.find((long)eff_orbit) != ZNArtag.end()) {
      for (auto zbc : ZNArtag[(long)eff_orbit]) {
        if ((int)(zbc / mStrobe) == (int)(eff_bc / mStrobe)) {
          isZNArtagged = true;
          Counters->Fill(4);
          break;
        }
      }
    }

    bool isZNAltagged = false;
    if (ZNAltag.find((long)eff_orbit) != ZNAltag.end()) {
      for (auto zbc : ZNAltag[(long)eff_orbit]) {
        if ((int)(zbc / mStrobe) == (int)(eff_bc / mStrobe)) {
          isZNAltagged = true;
          Counters->Fill(4);
          break;
        }
      }
    }

    bool isZNCrtagged = false;
    if (ZNCrtag.find((long)eff_orbit) != ZNCrtag.end()) {
      for (auto zbc : ZNCrtag[(long)eff_orbit]) {
        if ((int)(zbc / mStrobe) == (int)(eff_bc / mStrobe)) {
          isZNCrtagged = true;
          Counters->Fill(10);
          break;
        }
      }
    }

    bool isZNCltagged = false;
    if (ZNCltag.find((long)eff_orbit) != ZNCltag.end()) {
      for (auto zbc : ZNCltag[(long)eff_orbit]) {
        if ((int)(zbc / mStrobe) == (int)(eff_bc / mStrobe)) {
          isZNCltagged = true;
          Counters->Fill(10);
          break;
        }
      }
    }

    TZDCtag = 0;
    TZDCtag |= (isZNArtagged << 0);
    TZDCtag |= (isZNAltagged << 1);
    TZDCtag |= (isZNCrtagged << 2);
    TZDCtag |= (isZNCltagged << 3);

    // std::pair<int,int> closestbkg = findclosestbkg(Torbit, Tbc, ZDCAtag);

    // Tclosest_low = closestbkg.first;
    // Tclosest_up = closestbkg.second;

    // if (TZDCtag > 0) LOGP(info,"ZDCTag {}. Closest A-side bkg: asking {}/{}, returning (low,up) = ({},{})",TZDCtag,Torbit,Tbc,Tclosest_low,Tclosest_up);
    // else LOGP(info,"DEBUGcl bkg: asking {}/{}, returning {} , {}. ZDCTag {}",Torbit,Tbc,Tclosest_low,Tclosest_up,TZDCtag);

    if (TZDCtag > 0)
      LOGP(info, "ZDC tag with mask {}: ZNAright = {} - ZNAleft = {} - ZNCright = {} - ZNCleft = {}",
           TZDCtag, (TZDCtag >> 0) & 1, (TZDCtag >> 1) & 1, (TZDCtag >> 2) & 1, (TZDCtag >> 3) & 1);

    std::set<int> AvailableChips{};
    std::vector<std::set<int>> AvailableChipsInStave(12 + 16 + 20, std::set<int>{});
    std::map<int, std::vector<int>> MAPsize{}; // MAP[chip] = {list if sizes}
    std::map<int, std::vector<int>> MAPcols{}; // MAP[chip] = {list of column span}
    std::map<int, int> MAPntarget{};           // MAP[chip] = number of bad clusters in chip
    std::map<int, int> MAPcoo_mincol{};        // MAP[chip] = minimum of column coordinate
    std::map<int, int> MAPcoo_maxcol{};        // MAP[chip] = maximum of (column coordinate + colspan)

    int ntarget_in_rof = 0;
    for (int iclus = 0; iclus < clustersInRof.size(); iclus++) {

      const auto& compClus = clustersInRof[iclus];

      auto chipid = compClus.getSensorID();

      if (ChipToLayer(chipid) > 2) {
        continue;
      }

      ChipSeenInThisROF[chipid] = true;

      int coo_col = (int)compClus.getCol();
      int coo_row = (int)compClus.getRow();

      auto patti = patternsInRof[iclus];
      int npix = patti.getNPixels();
      int colspan = patti.getColumnSpan();
      int rowspan = patti.getRowSpan();

      if (npix > 1) {
        ClusterShape->Fill(colspan, rowspan);
        double pphi = ChipToPhi((int)chipid);
        if (pphi < 3.15 && pphi > 3.14 / 2)
          ClusterShapeTO->Fill(colspan, rowspan);
      }

      bool newchip = AvailableChips.insert(chipid).second;
      if (newchip) {
        MAPsize[chipid] = std::vector<int>{};
        MAPcols[chipid] = std::vector<int>{};
        MAPntarget[chipid] = 0;
        MAPcoo_mincol[chipid] = coo_col;
        MAPcoo_maxcol[chipid] = coo_col + colspan;
      }

      MAPsize[chipid].push_back(npix);
      MAPcols[chipid].push_back(colspan);
      if (colspan > 127 && rowspan < 30) {
        MAPntarget[chipid] += 1;
	ntarget_in_rof++;
      }
      MAPcoo_mincol[chipid] = TMath::Min(MAPcoo_mincol[chipid], coo_col);
      MAPcoo_maxcol[chipid] = TMath::Max(MAPcoo_maxcol[chipid], coo_col + colspan);

      AvailableChipsInStave[(int)(chipid / 9)].insert(chipid);

    } // end of loop over clusters in rof

    if (ntarget_in_rof == 0 && mTFn > mSkimmedOnlyAfterTF + 2){  // extra 2 to avoid edge effects?
      // do not need extra computations for this rof since it will not be saved in any case
      continue;
    }

    Ritstag = 0;
    RitstagTO = 0;
    Rnchip = RnchipTO = 0;
    RnchipafterTO = 0;
    Rnchipafter = 0;
    if (inTFROFcounter > 0) {
      for (int ic = 0; ic < N_CHIP_IB; ic++) {
        RnchipafterTO += (int)ChipSeenInLastROF[ic];
      }
    } else {
      RnchipafterTO = -1;
    }

    for (int ic : AvailableChips) {

      Rnchip++;
      Tchip = ic;

      if (inTFROFcounter < 1)
        Tmissingafter = -1;
      else if (ChipSeenInLastROF[ic])
        Tmissingafter = 0;
      else
        Tmissingafter = 1;

      if (inTFROFcounter < 2)
        Tmissingafter2 = -1;
      else if (ChipSeenInLast2ROF[ic])
        Tmissingafter2 = 0;
      else
        Tmissingafter2 = 1;

      Tphi = ChipToPhi(ic);

      if (Tphi < 3.15 && Tphi > 3.14 / 2)
        RnchipTO++;

      Tnclus = MAPsize[ic].size();
      Tmincol = MAPcoo_mincol[ic];
      Tmaxcol = MAPcoo_maxcol[ic];

      std::sort(MAPsize[ic].begin(), MAPsize[ic].end(), std::greater<int>());
      // std::sort(MAPcols[ic].begin(), MAPcols[ic].end(), std::greater<int>());

      Tnhit = Tnclus_s20 = Tnclus_s100 = Tnclus_s150 = 0;
      Tnhit1 = Tnhit10 = 0.;
      Tnclus_c20 = Tnclus_c100 = Tnclus_c128 = 0;
      Tnclus_target = MAPntarget[ic];
      Tnhit_no1pix = 0;
      Tnclus_no1pix = 0;

      int nhit_no1pix = 0;
      int n2hit = 0, n2hit_no1pix = 0;

      int nclus10 = 0, nclus1 = 0;

      for (int nh : MAPsize[ic]) {

        Tnhit += nh;

        n2hit += nh * nh;

        if (nh > 1) {
          Tnhit_no1pix += nh;
          Tnclus_no1pix += 1;
          n2hit_no1pix += nh * nh;
        }

        if (nclus10 < 10) {
          nclus10++;
          Tnhit10 += 1. * nh;
        }

        if (nclus1 < 1) {
          nclus1++;
          Tnhit1 += 1. * nh;
        }

        Tnclus_s20 += (nh >= 20);
        Tnclus_s100 += (nh >= 100);
        Tnclus_s150 += (nh >= 150);
      }

      Tnhit10 = (nclus10 == 0) ? 0. : 1. * Tnhit10 / nclus10;
      // computing STD as sqrt of E[X^2] - E[X]^2
      Tstdhit = (Tnclus == 0) ? 0. : TMath::Sqrt(1. * n2hit / Tnclus - TMath::Power(1. * Tnhit / Tnclus, 2));
      Tstdhit_no1pix = (Tnclus_no1pix == 0) ? 0. : TMath::Sqrt(1. * n2hit_no1pix / Tnclus_no1pix - TMath::Power(1. * Tnhit_no1pix / Tnclus_no1pix, 2));

      for (int nc : MAPcols[ic]) {
        Tnclus_c20 += (nc >= 20);
        Tnclus_c100 += (nc >= 100);
        Tnclus_c128 += (nc >= 128);
      }

      // counting events based on arbitrary definition of ITS-tagged
      if (Ritstag == 0 && Tnclus_target > 0) {
        Counters->Fill(5);
        if (TZDCtag > 0)
          Counters->Fill(7);
        Ritstag = 1;
      }
      if (RitstagTO == 0 && Tnclus_target > 0 && Tphi < 3.15 && Tphi > 3.14 / 2) {
        Counters->Fill(6);
        if (TZDCtag > 0)
          Counters->Fill(8);
        RitstagTO = 1;
      }

      // filling the Event Tree
      ITSChipEvtTree->Fill();
      if (Tnclus_target > 0){
	ITSChipEvtTargetTree->Fill();
      }

    } // end of loop over available chips

    for (int istave = 0; istave < (12 + 16 + 20); istave++) {

      if (AvailableChipsInStave[istave].size() == 0)
        continue;

      Snchip = 0;
      Sstave = istave;

      Snhit = Snclus = Snclus_s20 = Snclus_s100 = Snclus_s150 = 0;
      Snhit1 = Snhit10 = 0.;
      Snclus_c20 = Snclus_c100 = Snclus_c128 = Snclus_target = 0;

      std::vector<int> staveclustersizes{};
      std::vector<int> staveclusterscolumns{};
      Snclus_target = 0;

      for (int ic : AvailableChipsInStave[istave]) {

        staveclustersizes.insert(staveclustersizes.end(), MAPsize[ic].begin(), MAPsize[ic].end());
        staveclusterscolumns.insert(staveclusterscolumns.end(), MAPcols[ic].begin(), MAPcols[ic].end());
        Snclus_target += MAPntarget[ic];

        Snchip++;

        Sphi = ChipToPhi(ic);
      }

      Snclus = staveclustersizes.size();

      std::sort(staveclustersizes.begin(), staveclustersizes.end(), std::greater<int>());

      int nclus10 = 0, nclus1 = 0;

      for (int nh : staveclustersizes) {

        Snhit += nh;

        if (nclus10 < 10) {
          nclus10++;
          Snhit10 += 1. * nh;
        }

        if (nclus1 < 1) {
          nclus1++;
          Snhit1 += 1. * nh;
        }

        Snclus_s20 += (nh >= 20);
        Snclus_s100 += (nh >= 100);
        Snclus_s150 += (nh >= 150);
      }

      Snhit10 = (nclus10 == 0) ? 0. : 1. * Snhit10 / nclus10;

      /*
      for (int nc: staveclusterscolumns){
  Snclus_c20  += (nc >= 20);
  Snclus_c100 += (nc >= 100);
  Snclus_c128 += (nc >= 128);
      }
      */

      // filling the Event Tree
      ITSStaveEvtTree->Fill();

    } // end of loop over staves

    ITSROFEvtTree->Fill();

  } // end of loop over ROFs

  // plot interessanti:
  // Cluster size (averged or max) vs chip
  // Cluster length vs chip
  // Cluster prima d un evento in cui il chip è vuoto
  // numero di chip con cluster grande
  //    (x2 --> tutti gli eventi o con solo gli eventi ZDC-tagged)

  // Osservabile vs DeltaT (evento ITS e evento ZDC-tagged più vicino).

  // print tagged orbit-bc
  /*
  if (ZDCAtag.size()>0){
    for (auto ob : ZDCAtag){
      for (auto b : ob.second){
  LOGP(info,"ZNAC background on orbit/bc {}/{}",ob.first,b);
      }
    }
  }
  if (ZDCCtag.size()>0){
    for (auto ob : ZDCCtag){
      for (auto b : ob.second){
  LOGP(info,"ZNCC background on orbit/bc {}/{}",ob.first,b);
      }
    }
  }
  */
}

int ITSZDCAnomalyStudy::ChipToLayer(int chip)
{
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

double ITSZDCAnomalyStudy::ChipToPhi(int chip)
{

  int staveinlayer = (int)(chip / 9);
  for (int il = 0; il < ChipToLayer(chip); il++) {
    staveinlayer -= NStaves[il];
  }

  return 2. * TMath::Pi() * (0.5 + staveinlayer) / NStaves[ChipToLayer(chip)];
}

int ITSZDCAnomalyStudy::bcdistance(long orb1, int bc1, long orb2, int bc2)
{ // returns "2" - "1"
  if (orb1 == orb2)
    return bc2 - bc1;
  long dOrb = orb2 - orb1;
  int dOBC = static_cast<int>(dOrb * 3564);
  if (dOrb > 0)
    return dOBC - bc1 + bc2;
  else
    return 0 - ITSZDCAnomalyStudy::bcdistance(orb2, bc2, orb1, bc1);
}

std::pair<long, int> ITSZDCAnomalyStudy::shiftbc(long orbit, int bc, int shift)
{
  // TO BE IMPLEMENTED
  return std::make_pair(orbit, bc);
}

std::pair<int, int> ITSZDCAnomalyStudy::findclosestbkg(long orb, int bc, std::map<long, std::set<int>> ZDC)
{
  if (ZDC.size() == 0) {
    return std::make_pair(999999, -999999);
  }
  std::vector<int> distances{-999999};
  // std::vector<long> debug_orbit{-999999};
  // std::vector<int> debug_bc{-999999};
  for (auto ob : ZDC) {
    for (auto b : ob.second) {
      distances.push_back(bcdistance(orb, bc, ob.first, b));
      // debug_orbit.push_back(ob.first);
      // debug_bc.push_back(b);
    }
  }

  distances.push_back(999999);
  // debug_orbit.push_back(-999999);
  // debug_bc.push_back(-999999);

  int toret1, toret2;

  int i = 0;
  for (i = 0; i < distances.size() - 1; i++) {
    if (distances[i] < 0 && distances[i + 1] >= 0) {
      toret1 = distances[i];
      toret2 = distances[i + 1];
      break;
    }
  }

  // LOGP(info,"Closest bkg function. Asked {}/{} , returning distances {}({}/{}) , {}({}/{})", orb, bc, toret1, debug_orbit[i], debug_bc[i], toret2, debug_orbit[i+1], debug_bc[i+1]);

  if (toret1 == -999999)
    toret1 = 999999;
  if (toret2 == 999999)
    toret2 = -999999;

  return std::make_pair(toret1, toret2);
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
  // dataRequest->requestTracks(GTrackID::getSourcesMask("ZDC"), useMC); //working

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
