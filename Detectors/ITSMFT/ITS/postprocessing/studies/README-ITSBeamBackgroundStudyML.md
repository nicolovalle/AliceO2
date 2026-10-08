# ITS beam background study with unsupervised ML

`ITSBeamBackgroundStudyML` analyses inner-barrel layers 0–2 (chips 0–431).
It supports monolithic and per-layer/staggered ITS clusters. No labelled examples,
separate training file, or clean reference run is required.

The study collects a uniform reservoir of unlabelled clusters **across the whole
input stream**, independently for each layer. At end of stream it fits one
Isolation Forest per layer, then replays the stored features to score every
inner-barrel cluster, including early clusters. It is an offline study: final
scores and counts become available after end of stream. ROOT streams features
to disk; reservoir and forest sizes are bounded. Decoded patterns occupy memory
for one input slot of a timeframe.

The six features, in persisted order, are:

1. `log1p(npixels)`
2. `log1p(colspan)`
3. `log1p(rowspan)`
4. `npixels / (colspan * rowspan)` (bounding-box fill fraction)
5. `log(colspan / rowspan)`
6. `log1p(chipClusters)` (cluster count on that chip in that ROF)

Chip identity, detector coordinates, the legacy selection and external detector
tags are not model inputs. Tree training uses 80% of the reservoir; the remaining
20% calibrates the score threshold. Each tree samples without replacement from
the training partition. Larger scores indicate shorter average isolation paths,
following [Liu, Ting and Zhou, ICDM 2008](https://doi.org/10.1109/ICDM.2008.17).
Constant features are skipped when selecting a split; identical patterns remain
in a leaf with the expected-path-length correction.

The default selection is `score > threshold`, where the threshold is the 99th
percentile of held-out scores. This is a configurable tail selection policy,
not a measured false-positive rate or beam-background probability. Discrete/tied
scores may select fewer clusters than the requested tail fraction. Unsupervised
learning assumes anomalies are rare: common contamination, changes in operating
conditions and rare legitimate patterns can affect results. A feature that is
constant in the fitted sample cannot contribute a split for detecting a later
change. The legacy `colspan >= 128 && rowspan <= 29` flag is retained solely for
comparison. This task does not perform ZDC correlation or identify the physical
origin of an anomaly.

## Run

In an O2 environment, from a directory containing ITS cluster input files:

```sh
o2-its-standalone-postprocessing-workflow --its-beambkg-ml --disable-mc -b
```

The selection flag is `--its-beambkg-ml`; the internal DPL processor is named
`its-beambkg-ml-task`. DPL also registers processor names as options for grouped
settings, so the processor name must differ from the boolean selection flag.

Task options:

| Option | Default | Meaning |
| --- | --- | --- |
| `--ml-output` | `its_beambkg_ml.root` | ROOT feature, score and summary output |
| `--ml-reservoir-size` | `4096` | Maximum sampled clusters per layer; 64–65536 |
| `--ml-trees` | `100` | Trees per layer; 1–1024 |
| `--ml-sample-size` | `256` | Points per tree; 2–4096, capped to training partition size |
| `--ml-tail-fraction` | `0.01` | Calibration upper tail; strictly between 0 and 0.5 |
| `--ml-seed` | `42` | Nonnegative seed; layer number is added |
| `--ml-model-input` | empty | Optional previously saved model; skips fitting on this run |
| `--ml-model-output` | empty | Optional text model file for reuse |

At least 64 sampled clusters are needed to fit a layer. Insufficient layers are
explicitly marked unscored (`modelReady=false`, `score=threshold=-1`); they are
not counted as ML anomalies. An empty run still produces valid empty trees.
Fixed seed and fixed input order reproduce a model within the same C++ standard
library implementation. Model files are versioned and preserve feature order,
layer order, thresholds and tree populations. All three layer models are also
embedded in the ROOT output as `isolationForestModels`. The study uses the
topology dictionary supplied through DPL conditions for the input run.

### Comparing the legacy and ML selections

The boolean `clusters.legacyTarget` is true when the decoded cluster pattern has
`colspan >= 128 && rowspan <= 29`, matching the original beam-background study's
shape cuts. It is stored for comparison and is not used to train or score the
model. The ML selection is `scores.isAnomaly`, true for a scored cluster whose
`score > threshold`; `scores` is an entry-aligned friend of `clusters`.
The two selections can agree or disagree. Summing `legacyTarget` counts clusters
selected by the original cuts; summing `isAnomaly` counts the ML-selected clusters.

## ROOT output and counting

| Object | Contents |
| --- | --- |
| `clusters` | One entry per IB cluster: TF/ROF/slot, orbit/BC, chip/layer, cluster index, raw dimensions, occupancy, legacy flag, six features |
| `scores` | Entry-aligned friend of `clusters`: `modelReady`, `score`, `threshold`, `isAnomaly` |
| `rofs` | One entry per input ROF (per IB slot in staggered mode), including empty ROFs |
| `rofScores` | Entry-aligned friend of `rofs`: `nclus`, `nscored`, `nanomaly`, `nclus_target`, `maxScore` |
| `chipevtML` | One entry per nonempty chip–ROF pair, with identifiers and the same count branches |
| `Counters` | TFs, ROFs, IB clusters, scored clusters, ML anomalies, affected chip–ROF pairs, affected ROFs |
| `anomalyScores` | Score distributions per layer |
| `configuration`, `featureSchema` | Parameters and feature definition |

Sum `nanomaly` over **one** of `chipevtML` or `rofScores` to count anomalous
clusters. Count entries satisfying `nanomaly > 0` to count affected chip–ROF pairs
or ROFs respectively. `nclus_target` is the legacy shape-cut count, not the ML
count. ROF identity is `(TFprogress, slot, rof)`; staggered layer ROFs are
independent and are not collapsed into a common detector event.

The original `ITSBeamBackgroundStudy` and its workflow option remain available
for comparison. The new implementation requires no additional ML dependency.
