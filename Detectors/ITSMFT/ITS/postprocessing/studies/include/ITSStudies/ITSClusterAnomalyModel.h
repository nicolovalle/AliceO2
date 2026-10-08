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

#ifndef O2_ITS_CLUSTER_ANOMALY_MODEL_H
#define O2_ITS_CLUSTER_ANOMALY_MODEL_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <istream>
#include <limits>
#include <numeric>
#include <ostream>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace o2::its::study
{

// Feature order is part of the persisted model contract (version 1).
// Models are fitted separately for each IB layer; chip ID and the legacy cut are not features.
using ITSClusterFeatures = std::array<float, 6>;

inline ITSClusterFeatures makeITSClusterFeatures(int pixels, int columns, int rows, int chipClusters)
{
  if (pixels < 1 || columns < 1 || rows < 1 || pixels > columns * rows || chipClusters < 1) {
    throw std::invalid_argument("Invalid ITS cluster dimensions or occupancy");
  }
  return {std::log1p(float(pixels)), std::log1p(float(columns)), std::log1p(float(rows)),
          float(pixels) / (columns * rows), std::log(float(columns) / rows), std::log1p(float(chipClusters))};
}

// Bounded reservoir sampling followed by an Isolation Forest, using the path-length
// score of Liu, Ting and Zhou (ICDM 2008, doi:10.1109/ICDM.2008.17).
// The reservoir is unlabelled input data, not a required clean reference dataset.
class ITSClusterAnomalyModel
{
 public:
  static constexpr int MaxTrees = 1024;
  static constexpr int MaxSampleSize = 4096;
  static constexpr int MaxReservoirSize = 65536;

  explicit ITSClusterAnomalyModel(uint32_t seed = 42) : mRandom(seed) {}

  void observe(const ITSClusterFeatures& features, size_t capacity)
  {
    validate(features);
    if (ready() || capacity < 64 || capacity > MaxReservoirSize) {
      throw std::invalid_argument("Invalid reservoir capacity or already fitted model");
    }
    ++mSeen;
    if (mReservoir.size() < capacity) {
      mReservoir.push_back(features);
    } else {
      const auto index = std::uniform_int_distribution<uint64_t>(0, mSeen - 1)(mRandom);
      if (index < capacity) {
        mReservoir[index] = features;
      }
    }
  }

  // Hold out 20% of the reservoir for threshold calibration. No calibration points
  // participate in building trees. The tail fraction is a selection policy, not a probability.
  bool fit(int trees, int sampleSize, double tailFraction)
  {
    if (ready() || trees < 1 || trees > MaxTrees || sampleSize < 2 || sampleSize > MaxSampleSize ||
        !std::isfinite(tailFraction) || tailFraction <= 0. || tailFraction >= 0.5) {
      throw std::invalid_argument("Invalid Isolation Forest configuration");
    }
    if (mReservoir.size() < 64) {
      return false;
    }
    mReferenceSize = mReservoir.size();
    std::shuffle(mReservoir.begin(), mReservoir.end(), mRandom);
    const int calibrationSize = mReferenceSize / 5;
    const int trainingSize = mReferenceSize - calibrationSize;
    mSampleSize = std::min(sampleSize, trainingSize);
    mPathCorrection.resize(mSampleSize + 1, 0.);
    double harmonic = 0.;
    for (int n = 2; n <= mSampleSize; ++n) {
      harmonic += 1. / (n - 1);
      mPathCorrection[n] = 2. * harmonic - 2. * (n - 1) / n;
    }
    const int maxDepth = int(std::ceil(std::log2(mSampleSize)));
    std::vector<int> indices(trainingSize);
    std::iota(indices.begin(), indices.end(), 0);
    mForest.resize(trees);
    for (auto& tree : mForest) {
      std::shuffle(indices.begin(), indices.end(), mRandom);
      std::vector<int> sample(indices.begin(), indices.begin() + mSampleSize);
      tree.reserve(2 * mSampleSize - 1);
      build(tree, sample, 0, maxDepth);
    }
    std::vector<double> scores;
    scores.reserve(calibrationSize);
    for (size_t i = trainingSize; i < mReferenceSize; ++i) {
      scores.push_back(score(mReservoir[i]));
    }
    std::sort(scores.begin(), scores.end());
    const auto quantile = std::min(scores.size() - 1, size_t(std::ceil((1. - tailFraction) * scores.size()) - 1));
    mThreshold = scores[quantile];
    mReservoir.clear();
    mReservoir.shrink_to_fit();
    return true;
  }

  bool ready() const { return !mForest.empty(); }
  double threshold() const { return mThreshold; }
  uint64_t seen() const { return mSeen; }
  size_t referenceSize() const { return mReferenceSize; }
  size_t reservoirSize() const { return mReservoir.size(); }
  int sampleSize() const { return mSampleSize; }
  size_t trees() const { return mForest.size(); }

  double score(const ITSClusterFeatures& features) const
  {
    validate(features);
    if (!ready()) {
      throw std::logic_error("Scoring an unfitted ITS anomaly model");
    }
    double path = 0.;
    for (const auto& tree : mForest) {
      int index = 0, depth = 0;
      while (tree[index].feature >= 0) {
        const auto& node = tree[index];
        index = features[node.feature] < node.split ? node.left : node.right;
        ++depth;
      }
      path += depth + mPathCorrection[tree[index].size];
    }
    return std::exp2(-path / (mForest.size() * mPathCorrection[mSampleSize]));
  }

  void write(std::ostream& out) const
  {
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << ready() << ' ' << mSeen << ' ' << mReferenceSize << ' ' << mThreshold << ' ' << mSampleSize << ' ' << mForest.size() << '\n';
    for (const auto& tree : mForest) {
      out << tree.size() << '\n';
      // Nodes are stored in preorder; child indices are reconstructed on load.
      for (const auto& node : tree) {
        out << node.feature << ' ' << node.split << ' ' << node.size << '\n';
      }
    }
    if (!out) {
      throw std::runtime_error("Cannot write ITS anomaly model");
    }
  }

  void read(std::istream& in)
  {
    ITSClusterAnomalyModel loaded;
    int fitted;
    size_t nTrees;
    if (!(in >> fitted >> loaded.mSeen >> loaded.mReferenceSize >> loaded.mThreshold >> loaded.mSampleSize >> nTrees) ||
        (fitted != 0 && fitted != 1) || !std::isfinite(loaded.mThreshold) || loaded.mThreshold < 0. || loaded.mThreshold > 1. ||
        loaded.mReferenceSize > MaxReservoirSize || loaded.mReferenceSize > loaded.mSeen ||
        (fitted && (nTrees < 1 || nTrees > MaxTrees || loaded.mSampleSize < 2 || loaded.mSampleSize > MaxSampleSize || loaded.mSampleSize > int(loaded.mReferenceSize))) ||
        (!fitted && (nTrees != 0 || loaded.mSampleSize != 0))) {
      throw std::runtime_error("Invalid ITS anomaly model header");
    }
    loaded.mPathCorrection.resize(loaded.mSampleSize + 1, 0.);
    double harmonic = 0.;
    for (int n = 2; n <= loaded.mSampleSize; ++n) {
      harmonic += 1. / (n - 1);
      loaded.mPathCorrection[n] = 2. * harmonic - 2. * (n - 1) / n;
    }
    loaded.mForest.resize(nTrees);
    for (auto& tree : loaded.mForest) {
      int nodes;
      if (!(in >> nodes) || nodes < 1 || nodes > 2 * loaded.mSampleSize - 1) {
        throw std::runtime_error("Invalid isolation tree size");
      }
      tree.reserve(nodes);
      loaded.readNode(in, tree, nodes, 0);
      if (tree.size() != size_t(nodes) || tree.front().size != loaded.mSampleSize) {
        throw std::runtime_error("Incomplete isolation tree");
      }
    }
    *this = std::move(loaded);
  }

 private:
  struct Node {
    int feature = -1;
    float split = 0.;
    int size = 0;
    int left = -1, right = -1;
  };

  static void validate(const ITSClusterFeatures& features)
  {
    if (!std::all_of(features.begin(), features.end(), [](float x) { return std::isfinite(x); })) {
      throw std::invalid_argument("Nonfinite ITS cluster features");
    }
  }

  int build(std::vector<Node>& tree, const std::vector<int>& sample, int depth, int maxDepth)
  {
    const int index = tree.size();
    tree.push_back(Node{-1, 0., int(sample.size()), -1, -1});
    if (sample.size() < 2 || depth >= maxDepth) {
      return index;
    }
    ITSClusterFeatures low = mReservoir[sample.front()], high = low;
    for (int i : sample) {
      for (size_t f = 0; f < low.size(); ++f) {
        low[f] = std::min(low[f], mReservoir[i][f]);
        high[f] = std::max(high[f], mReservoir[i][f]);
      }
    }
    std::vector<int> varying;
    for (size_t f = 0; f < low.size(); ++f) {
      if (low[f] < high[f]) {
        varying.push_back(f);
      }
    }
    if (varying.empty()) {
      return index;
    }
    const int feature = varying[std::uniform_int_distribution<size_t>(0, varying.size() - 1)(mRandom)];
    // Keep the split strictly above min and at most max, including adjacent float values.
    const float split = std::max(std::nextafter(low[feature], high[feature]), std::uniform_real_distribution<float>(low[feature], high[feature])(mRandom));
    std::vector<int> left, right;
    for (int i : sample) {
      (mReservoir[i][feature] < split ? left : right).push_back(i);
    }
    tree[index].feature = feature;
    tree[index].split = split;
    const int leftIndex = build(tree, left, depth + 1, maxDepth);
    const int rightIndex = build(tree, right, depth + 1, maxDepth);
    tree[index].left = leftIndex;
    tree[index].right = rightIndex;
    return index;
  }

  int readNode(std::istream& in, std::vector<Node>& tree, int nodes, int depth)
  {
    Node node;
    if (tree.size() >= size_t(nodes) || depth > int(std::ceil(std::log2(mSampleSize))) ||
        !(in >> node.feature >> node.split >> node.size) || node.feature < -1 || node.feature >= int(ITSClusterFeatures{}.size()) ||
        !std::isfinite(node.split) || node.size < 1 || node.size > mSampleSize) {
      throw std::runtime_error("Invalid isolation tree node");
    }
    const int index = tree.size();
    tree.push_back(node);
    if (node.feature >= 0) {
      const int left = readNode(in, tree, nodes, depth + 1);
      const int right = readNode(in, tree, nodes, depth + 1);
      if (tree[left].size + tree[right].size != node.size) {
        throw std::runtime_error("Inconsistent isolation tree populations");
      }
      tree[index].left = left;
      tree[index].right = right;
    }
    return index;
  }

  std::mt19937 mRandom;
  std::vector<ITSClusterFeatures> mReservoir;
  uint64_t mSeen = 0;
  size_t mReferenceSize = 0;
  int mSampleSize = 0;
  double mThreshold = 1.;
  std::vector<double> mPathCorrection;
  std::vector<std::vector<Node>> mForest;
};
} // namespace o2::its::study

#endif
