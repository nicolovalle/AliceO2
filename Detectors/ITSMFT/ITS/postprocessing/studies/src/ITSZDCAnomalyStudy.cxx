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
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  std::shared_ptr<DataRequest> mDataRequest;
  bool mUseMC;
  size_t mTFn = 0;
  const o2::itsmft::TopologyDictionary *mDict = nullptr;
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
  mgr.setTimestamp(o2::ccdb::getCurrentTimestamp());
  mDict = mgr.get<o2::itsmft::TopologyDictionary>("ITS/Calib/ClusterDictionary");
  
}

void ITSZDCAnomalyStudy::endOfStream(EndOfStreamContext&)
{
  LOGP(info, "End of stream for ITSZDCAnomalyStudy");
}

void ITSZDCAnomalyStudy::run(ProcessingContext& pc)
{
  LOGP(info, "Running ITSZDCAnomalyStudy on TF: {}", mTFn++);
  o2::globaltracking::RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest.get());
  // updateTimeDependentParams(pc);
  process(recoData);
}

void ITSZDCAnomalyStudy::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  //   if (matcher == ConcreteDataMatcher("ITS", "CLUSDICT", 0)) {
  //     setClusterDictionary((const o2::itsmft::TopologyDictionary*)obj);
  //     return;
  //   }
}

// Custom area
void ITSZDCAnomalyStudy::process(o2::globaltracking::RecoContainer& recoData)
{
  o2::zdc::RecEventFlat ev;

  LOGP(info, "Processing RecoContainer");
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
LOGP(info, "sizeof ITS RC : {}, {}, {}", clusArr.size(), clusPatt.size(), rofRecVec.size());

  ev.init(RecBC, Energy, TDCData, Info2);
  while (ev.next()) {
    // cout<<(int)(ev.getNTDC())<<endl;
    int32_t itdc = 0;
    int nhit = ev.NtdcV(itdc);
    for (int32_t ipos = 0; ipos < nhit; ipos++) {
      double mytdc = o2::zdc::FTDCVal * ev.TDCVal[itdc][ipos];
      // ht[ich]->Fill(mytdc);
      if (itdc == o2::zdc::TDCZNAC && mytdc > 5.7 && mytdc < 8.7) {
        LOGF(info, "ZNAC background %f hit on %u.%04u\n", mytdc, ev.ir.orbit, ev.ir.bc);
      }
    }
  }
}

// getter
DataProcessorSpec getITSZDCAnomalyStudy(mask_t srcClustersMask, bool useMC)
{
  std::vector<OutputSpec> outputs;
  auto dataRequest = std::make_shared<DataRequest>();
  dataRequest->requestClusters(srcClustersMask, useMC);
  dataRequest->requestTracks(GTrackID::getSourcesMask("ZDC"), useMC);

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
