// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "ITSStudies/ITSBeamBackgroundStudyML.h"
#include "ITSStudies/ITSClusterAnomalyModel.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DataFormatsITSMFT/ClusterPattern.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "Framework/Task.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/Logger.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TObjString.h>
#include <TTree.h>

#include <array>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace o2::its::study
{
using namespace o2::framework;
using namespace o2::globaltracking;

class ITSBeamBackgroundStudyML : public Task
{
 public:
  explicit ITSBeamBackgroundStudyML(std::shared_ptr<DataRequest> request) : mDataRequest(std::move(request)) {}
  void init(InitContext&) final;
  void run(ProcessingContext&) final;
  void endOfStream(EndOfStreamContext&) final;

 private:
  static constexpr int NLayer = 3;
  static constexpr int NChip = o2::itsmft::ChipMappingITS::getFirstChipsOnLayer(NLayer);
  static constexpr const char* ModelHeader = "ITSBeamBackgroundStudyML 1 3";

  struct ROFRow {
    uint64_t tf = 0;
    uint32_t orbit = 0;
    int bc = 0, slot = 0, rof = 0, clusters = 0;
    Long64_t firstCluster = 0;
  } mROF;

  struct ClusterRow {
    int chip = 0, layer = 0, index = 0, row = 0, col = 0;
    int pixels = 0, columns = 0, rows = 0, chipClusters = 0;
    bool legacyTarget = false;
    ITSClusterFeatures features{};
  } mCluster;

  struct ScoreRow {
    bool ready = false, anomaly = false;
    double score = -1., threshold = -1.;
  } mScore;

  struct Counts {
    int clusters = 0, scored = 0, anomalies = 0, legacy = 0;
    double maxScore = -1.;
  } mChipCounts, mROFCounts;

  void process(RecoContainer&, const o2::itsmft::TopologyDictionary&);
  void scoreClusters();
  void writeModels(std::ostream&) const;
  static void fill(TTree* tree);
  TTree* makeTree(const char* name, const char* title);
  void addROFBranches(TTree*);
  static void addCountBranches(TTree*, Counts&);

  std::shared_ptr<DataRequest> mDataRequest;
  std::array<ITSClusterAnomalyModel, NLayer> mModels;
  std::unique_ptr<TFile> mOutput;
  TTree *mClusters = nullptr, *mROFs = nullptr;
  TH1D* mCounters = nullptr;
  TH2D* mScoreHistogram = nullptr;
  int mReservoirSize = 4096, mTrees = 100, mSampleSize = 256;
  double mTailFraction = 0.01;
  std::string mModelInput, mModelOutput;
  uint64_t mTFCount = 0;
};

TTree* ITSBeamBackgroundStudyML::makeTree(const char* name, const char* title)
{
  auto* tree = new TTree(name, title);
  tree->SetDirectory(mOutput.get());
  tree->SetAutoFlush(-8 * 1024 * 1024);
  tree->SetAutoSave(-64 * 1024 * 1024);
  return tree;
}

void ITSBeamBackgroundStudyML::fill(TTree* tree)
{
  if (tree->Fill() < 0) {
    throw std::runtime_error("Cannot write ITS ML study tree");
  }
}

void ITSBeamBackgroundStudyML::addROFBranches(TTree* tree)
{
  tree->Branch("TFprogress", &mROF.tf, "TFprogress/l");
  tree->Branch("orbit", &mROF.orbit, "orbit/i");
  tree->Branch("bc", &mROF.bc, "bc/I");
  tree->Branch("slot", &mROF.slot, "slot/I");
  tree->Branch("rof", &mROF.rof, "rof/I");
}

void ITSBeamBackgroundStudyML::addCountBranches(TTree* tree, Counts& counts)
{
  tree->Branch("nclus", &counts.clusters, "nclus/I");
  tree->Branch("nscored", &counts.scored, "nscored/I");
  tree->Branch("nanomaly", &counts.anomalies, "nanomaly/I");
  tree->Branch("nclus_target", &counts.legacy, "nclus_target/I");
  tree->Branch("maxScore", &counts.maxScore, "maxScore/D");
}

void ITSBeamBackgroundStudyML::init(InitContext& ic)
{
  mReservoirSize = ic.options().get<int>("ml-reservoir-size");
  mTrees = ic.options().get<int>("ml-trees");
  mSampleSize = ic.options().get<int>("ml-sample-size");
  mTailFraction = ic.options().get<double>("ml-tail-fraction");
  const int seed = ic.options().get<int>("ml-seed");
  mModelInput = ic.options().get<std::string>("ml-model-input");
  mModelOutput = ic.options().get<std::string>("ml-model-output");
  const auto output = ic.options().get<std::string>("ml-output");
  if (mReservoirSize < 64 || mReservoirSize > ITSClusterAnomalyModel::MaxReservoirSize ||
      mTrees < 1 || mTrees > ITSClusterAnomalyModel::MaxTrees ||
      mSampleSize < 2 || mSampleSize > ITSClusterAnomalyModel::MaxSampleSize ||
      !std::isfinite(mTailFraction) || mTailFraction <= 0. || mTailFraction >= 0.5 || seed < 0 || output.empty() ||
      output == mModelInput || output == mModelOutput || (!mModelInput.empty() && mModelInput == mModelOutput)) {
    throw std::invalid_argument("Invalid ITSBeamBackgroundStudyML configuration or overlapping output paths");
  }
  for (int layer = 0; layer < NLayer; ++layer) {
    mModels[layer] = ITSClusterAnomalyModel(uint32_t(seed) + layer);
  }
  if (!mModelInput.empty()) {
    std::ifstream in(mModelInput);
    std::string header;
    if (!std::getline(in, header) || header != ModelHeader) {
      throw std::runtime_error("ITS ML model has an incompatible feature schema or layer count");
    }
    for (auto& model : mModels) {
      model.read(in);
    }
    in >> std::ws;
    if (!in.eof()) {
      throw std::runtime_error("Unexpected trailing data in ITS ML model");
    }
    LOGP(info, "Loaded ITS inner-barrel anomaly models from {}", mModelInput);
  }

  mOutput.reset(TFile::Open(output.c_str(), "RECREATE"));
  if (!mOutput || mOutput->IsZombie()) {
    throw std::runtime_error("Cannot open ITS ML study ROOT output");
  }
  mOutput->cd();
  mClusters = makeTree("clusters", "Inner-barrel cluster features; scores are in the scores friend");
  addROFBranches(mClusters);
  mClusters->Branch("chip", &mCluster.chip, "chip/I");
  mClusters->Branch("layer", &mCluster.layer, "layer/I");
  mClusters->Branch("clusterIndex", &mCluster.index, "clusterIndex/I");
  mClusters->Branch("row", &mCluster.row, "row/I");
  mClusters->Branch("col", &mCluster.col, "col/I");
  mClusters->Branch("npixels", &mCluster.pixels, "npixels/I");
  mClusters->Branch("colspan", &mCluster.columns, "colspan/I");
  mClusters->Branch("rowspan", &mCluster.rows, "rowspan/I");
  mClusters->Branch("chipClusters", &mCluster.chipClusters, "chipClusters/I");
  mClusters->Branch("legacyTarget", &mCluster.legacyTarget, "legacyTarget/O");
  mClusters->Branch("features", mCluster.features.data(), "features[6]/F");
  mROFs = makeTree("rofs", "Input ROFs, including ROFs with no inner-barrel clusters");
  addROFBranches(mROFs);
  mROFs->Branch("firstCluster", &mROF.firstCluster, "firstCluster/L");
  mROFs->Branch("nclusIB", &mROF.clusters, "nclusIB/I");

  mCounters = new TH1D("Counters", "ITS ML study counts", 7, 0.5, 7.5);
  const std::array<const char*, 7> labels{"TF", "ROF", "IB clusters", "Scored clusters", "ML anomalies", "Affected chip-ROFs", "Affected ROFs"};
  for (int i = 0; i < int(labels.size()); ++i) {
    mCounters->GetXaxis()->SetBinLabel(i + 1, labels[i]);
  }
  mScoreHistogram = new TH2D("anomalyScores", "Isolation Forest score;Layer;Score", NLayer, -0.5, NLayer - 0.5, 100, 0., 1.);
  TNamed("featureSchema", "v1: log1p(npixels), log1p(colspan), log1p(rowspan), fillFraction, log(colspan/rowspan), log1p(chipClusters)").Write();
  std::ostringstream config;
  config << "reservoir=" << mReservoirSize << ";trees=" << mTrees << ";sample=" << mSampleSize
         << ";tailFraction=" << mTailFraction << ";seed=" << seed << ";modelInput=" << mModelInput;
  TNamed("configuration", config.str().c_str()).Write();
  LOGP(info, "ITS ML study: collecting unlabelled IB clusters; fitting and scoring at end of stream");
}

void ITSBeamBackgroundStudyML::run(ProcessingContext& pc)
{
  RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest);
  // Use the condition supplied for this run, rather than a dictionary at the current wall-clock time.
  const auto dictionary = pc.inputs().get<o2::itsmft::TopologyDictionary*>("cldictITS");
  if (!dictionary) {
    throw std::runtime_error("Missing ITS cluster topology dictionary");
  }
  ++mTFCount;
  mCounters->Fill(1);
  process(recoData, *dictionary);
}

void ITSBeamBackgroundStudyML::process(RecoContainer& recoData, const o2::itsmft::TopologyDictionary& dictionary)
{
  // Even OB explicit patterns must be consumed in monolithic input to preserve byte-stream alignment.
  const int slots = recoData.getITSPerLayer() ? NLayer : 1;
  for (int slot = 0; slot < slots; ++slot) {
    const auto clusters = recoData.getITSClusters(slot);
    const auto bytes = recoData.getITSClustersPatterns(slot);
    const auto rofs = recoData.getITSClustersROFRecords(slot);
    auto patternIt = bytes.begin();
    std::vector<o2::itsmft::ClusterPattern> patterns;
    patterns.reserve(clusters.size());
    for (const auto& cluster : clusters) {
      const auto id = cluster.getPatternID();
      o2::itsmft::ClusterPattern pattern;
      if (id != o2::itsmft::CompCluster::InvalidPatternID && id >= dictionary.getSize()) {
        throw std::runtime_error("ITS cluster pattern ID outside the topology dictionary");
      }
      if (id == o2::itsmft::CompCluster::InvalidPatternID || dictionary.isGroup(id)) {
        const auto remaining = bytes.end() - patternIt;
        if (remaining < 2) {
          throw std::runtime_error("Truncated ITS explicit cluster pattern header");
        }
        const int rows = patternIt[0], cols = patternIt[1];
        const int payload = (rows * cols + 7) / 8;
        if (rows < 1 || cols < 1 || rows > o2::itsmft::ClusterPattern::MaxRowSpan || cols > o2::itsmft::ClusterPattern::MaxColSpan || remaining < 2 + payload) {
          throw std::runtime_error("Invalid or truncated ITS explicit cluster pattern");
        }
        pattern.acquirePattern(patternIt);
      } else {
        pattern = dictionary.getPattern(id);
      }
      patterns.push_back(pattern);
    }
    if (patternIt != bytes.end()) {
      throw std::runtime_error("Unused bytes in ITS cluster pattern stream");
    }

    int rofIndex = 0;
    for (const auto& rof : rofs) {
      if (rof.getFirstEntry() < 0 || rof.getNEntries() < 0 ||
          size_t(rof.getFirstEntry()) > clusters.size() || size_t(rof.getNEntries()) > clusters.size() - rof.getFirstEntry()) {
        throw std::runtime_error("ITS ROF cluster range outside the cluster array");
      }
      mROF = {mTFCount, rof.getBCData().orbit, rof.getBCData().bc, slot, rofIndex++, 0, mClusters->GetEntries()};
      const auto rofClusters = rof.getROFData(clusters);
      const auto rofPatterns = rof.getROFData(patterns);
      std::array<int, NChip> occupancy{};
      for (const auto& cluster : rofClusters) {
        const int chip = cluster.getSensorID();
        if (chip >= o2::itsmft::ChipMappingITS::getNChips()) {
          throw std::runtime_error("Invalid ITS chip ID");
        }
        if (chip < NChip) {
          if (recoData.getITSPerLayer() && o2::itsmft::ChipMappingITS::getLayer(chip) != slot) {
            throw std::runtime_error("ITS chip layer does not match staggered input slot");
          }
          ++occupancy[chip];
        }
      }
      for (int i = 0; i < int(rofClusters.size()); ++i) {
        const auto& cluster = rofClusters[i];
        const int chip = cluster.getSensorID();
        if (chip >= NChip) {
          continue;
        }
        const auto& pattern = rofPatterns[i];
        mCluster.chip = chip;
        mCluster.layer = o2::itsmft::ChipMappingITS::getLayer(chip);
        mCluster.index = rof.getFirstEntry() + i;
        mCluster.row = cluster.getRow();
        mCluster.col = cluster.getCol();
        mCluster.pixels = pattern.getNPixels();
        mCluster.columns = pattern.getColumnSpan();
        mCluster.rows = pattern.getRowSpan();
        mCluster.chipClusters = occupancy[chip];
        mCluster.legacyTarget = mCluster.columns >= 128 && mCluster.rows <= 29;
        mCluster.features = makeITSClusterFeatures(mCluster.pixels, mCluster.columns, mCluster.rows, mCluster.chipClusters);
        if (mModelInput.empty()) {
          mModels[mCluster.layer].observe(mCluster.features, mReservoirSize);
        }
        fill(mClusters);
        ++mROF.clusters;
      }
      fill(mROFs);
      mCounters->Fill(2);
      mCounters->Fill(3, mROF.clusters);
    }
  }
  LOGP(info, "ITS ML study TF {}: {} IB cluster feature entries collected", mTFCount, mClusters->GetEntries());
}

void ITSBeamBackgroundStudyML::scoreClusters()
{
  // Flush feature baskets to disk before replay; memory use is bounded by ROOT baskets,
  // the per-layer reservoirs, and the forests, rather than the length of the run.
  mClusters->FlushBaskets();
  mROFs->FlushBaskets();
  mOutput->cd();
  auto* scores = makeTree("scores", "Entry-aligned cluster scores; -1 means layer model unavailable");
  scores->Branch("modelReady", &mScore.ready, "modelReady/O");
  scores->Branch("score", &mScore.score, "score/D");
  scores->Branch("threshold", &mScore.threshold, "threshold/D");
  scores->Branch("isAnomaly", &mScore.anomaly, "isAnomaly/O");
  auto* chipEvents = makeTree("chipevtML", "One entry per nonempty IB chip and ROF");
  addROFBranches(chipEvents);
  chipEvents->Branch("chip", &mCluster.chip, "chip/I");
  chipEvents->Branch("layer", &mCluster.layer, "layer/I");
  addCountBranches(chipEvents, mChipCounts);
  auto* rofScores = makeTree("rofScores", "Entry-aligned ROF counts; includes empty ROFs");
  addCountBranches(rofScores, mROFCounts);

  for (Long64_t iROF = 0; iROF < mROFs->GetEntries(); ++iROF) {
    if (mROFs->GetEntry(iROF) <= 0) {
      throw std::runtime_error("Cannot replay ITS ROF features");
    }
    std::array<Counts, NChip> counts{};
    mROFCounts = {};
    for (int i = 0; i < mROF.clusters; ++i) {
      if (mClusters->GetEntry(mROF.firstCluster + i) <= 0) {
        throw std::runtime_error("Cannot replay ITS cluster features");
      }
      const auto& model = mModels[mCluster.layer];
      mScore = {};
      mScore.ready = model.ready();
      if (mScore.ready) {
        mScore.score = model.score(mCluster.features);
        mScore.threshold = model.threshold();
        // Strict comparison prevents tied/discrete scores from forcing a fixed anomaly count.
        mScore.anomaly = mScore.score > mScore.threshold;
        mCounters->Fill(4);
        if (mScore.anomaly) {
          mCounters->Fill(5);
        }
        mScoreHistogram->Fill(mCluster.layer, mScore.score);
      }
      fill(scores);
      auto& chip = counts[mCluster.chip];
      ++chip.clusters;
      chip.scored += mScore.ready;
      chip.anomalies += mScore.anomaly;
      chip.legacy += mCluster.legacyTarget;
      chip.maxScore = std::max(chip.maxScore, mScore.score);
    }
    for (int chip = 0; chip < NChip; ++chip) {
      if (!counts[chip].clusters) {
        continue;
      }
      mCluster.chip = chip;
      mCluster.layer = o2::itsmft::ChipMappingITS::getLayer(chip);
      mChipCounts = counts[chip];
      fill(chipEvents);
      mROFCounts.clusters += mChipCounts.clusters;
      mROFCounts.scored += mChipCounts.scored;
      mROFCounts.anomalies += mChipCounts.anomalies;
      mROFCounts.legacy += mChipCounts.legacy;
      mROFCounts.maxScore = std::max(mROFCounts.maxScore, mChipCounts.maxScore);
      if (mChipCounts.anomalies) {
        mCounters->Fill(6);
      }
    }
    fill(rofScores);
    if (mROFCounts.anomalies) {
      mCounters->Fill(7);
    }
  }
  if (scores->GetEntries() != mClusters->GetEntries() || rofScores->GetEntries() != mROFs->GetEntries()) {
    throw std::runtime_error("ITS ML feature and score trees are misaligned");
  }
  mClusters->AddFriend(scores);
  mROFs->AddFriend(rofScores);
}

void ITSBeamBackgroundStudyML::writeModels(std::ostream& out) const
{
  out << ModelHeader << '\n';
  for (const auto& model : mModels) {
    model.write(out);
  }
}

void ITSBeamBackgroundStudyML::endOfStream(EndOfStreamContext&)
{
  for (int layer = 0; layer < NLayer; ++layer) {
    auto& model = mModels[layer];
    if (mModelInput.empty()) {
      model.fit(mTrees, mSampleSize, mTailFraction);
    }
    if (model.ready()) {
      LOGP(info, "ITS ML layer {}: {} input clusters, {} reservoir points, {} trees, sample size {}, threshold {}",
           layer, model.seen(), model.referenceSize(), model.trees(), model.sampleSize(), model.threshold());
    } else {
      LOGP(warning, "ITS ML layer {} has no fitted model (at least 64 clusters are needed); its clusters will be marked unscored", layer);
    }
  }
  scoreClusters();
  mOutput->cd();
  std::ostringstream models;
  writeModels(models);
  TObjString(models.str().c_str()).Write("isolationForestModels");
  if (!mModelOutput.empty()) {
    std::ofstream out(mModelOutput);
    writeModels(out);
    out.close();
    if (!out) {
      throw std::runtime_error("Cannot save ITS ML model file");
    }
  }
  if (mOutput->Write() < 0 || mOutput->TestBit(TFile::kWriteError)) {
    throw std::runtime_error("Cannot finalize ITS ML study ROOT output");
  }
  LOGP(info, "ITS ML study: {} clusters scored, {} anomalies; output {}", mCounters->GetBinContent(4), mCounters->GetBinContent(5), mOutput->GetName());
  mOutput->Close();
}

DataProcessorSpec getITSBeamBackgroundStudyML(bool useMC, bool itsStag)
{
  auto request = std::make_shared<DataRequest>();
  request->setITSPerLayer(itsStag);
  request->requestITSClusters(useMC);
  return DataProcessorSpec{
    // DPL registers a grouped-options argument using this name. It must differ
    // from the boolean workflow flag --its-beambkg-ml.
    "its-beambkg-ml-task",
    request->inputs,
    {},
    AlgorithmSpec{adaptFromTask<ITSBeamBackgroundStudyML>(request)},
    Options{
      {"ml-output", VariantType::String, "its_beambkg_ml.root", {"ROOT output for features, scores and counts"}},
      {"ml-reservoir-size", VariantType::Int, 4096, {"Maximum unlabelled reservoir points per IB layer (64..65536)"}},
      {"ml-trees", VariantType::Int, 100, {"Isolation trees per layer (1..1024)"}},
      {"ml-sample-size", VariantType::Int, 256, {"Training points per tree (2..4096; capped to available training data)"}},
      {"ml-tail-fraction", VariantType::Double, 0.01, {"Upper calibration score tail used for selection (0 < value < 0.5)"}},
      {"ml-seed", VariantType::Int, 42, {"Nonnegative random seed; layer number is added to it"}},
      {"ml-model-input", VariantType::String, "", {"Load a saved model instead of fitting on this run"}},
      {"ml-model-output", VariantType::String, "", {"Optionally save models for use on another run"}}}};
}
} // namespace o2::its::study
