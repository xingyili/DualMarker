# DualMarker

DualMarker: A multi-source fusion identification method for prognostic biomarkers of breast cancer based on dual-layer heterogeneous network.

## Description

![img](Structure.png)

## Getting Started

### Introduction

DualMarker is a biomarker recognition algorithm.

### Special dependencies

* **MATLAB** R2020a or later (R2021b+ recommended)
* **Statistics and Machine Learning Toolbox** (for standardization & classifier; you can also use `TreeBagger`)
* **Parallel Computing Toolbox** *(optional)* for faster repeated CV

> The repository includes implementations of `classRF_train/classRF_predict` and `AUC`. If you prefer MATLAB’s native classifier, replace them with `TreeBagger/predict` (an example is provided below).

---

### Overview

The repository is organized as follows:

* `DualMarkerMainCode.m`

  Main pipeline of  **DualMarker** :

  1. Build the heterogeneous network and compute RandNE embeddings;
  2. Derive prior scores from DisGeNET (gene–disease) and expression (t-test), normalize and fuse;
  3. Perform random walk / network propagation on the **dual-layer** network;
  4. Output final **biomarker ranking** in variable `geneRANK`.
* `KFoldCrossValidation.m`

  **Repeated 5-fold cross-validation** (default  **50 repeats × 5 folds** ) on **Top-100/150/200** features selected from `geneRANK`. Trains a classifier (random forest by default) and reports AUC for each repeat as well as the final mean AUCs.
* Core utilities

  `A_RWRplus.m`, `Network_Enhancement.m`, `getNormalizedMatrix_Heter.m`,

  `getRandNEemb_in.m` / `RandNE_*`, `Obtain_StatisticScore.m`, `TransitionFields.m`, `GS.m`, …
* `Data/`

  **One packaged `.mat` file** with everything needed to run DualMarker.

---

### Data

We provide a packaged file (example name):

```
Data/initavalue2.mat
```

This `.mat` file bundles all inputs needed to run DualMarker:

* **Heterogeneous network adjacency:** a single sparse block matrix encoding edges among genes, proteins, diseases, and GO terms (e.g., gene–gene, gene–protein, protein–protein, gene–disease, gene–GO). Used to build the dual-layer network and for propagation.
* **Expression matrix & gene list:** normalized gene expression values (samples × genes) with matching gene identifiers. Used to compute differential signals and for 5-fold cross-validation.
* **DisGeNET prior:** disease–gene relevance scores that act as prior seeds to guide the propagation and ranking.

### Running DualMarker

**1) Produce the biomarker ranking**

```matlab
addpath(genpath(pwd));                % add repo to MATLAB path
load('Data/DualMarker_demo.mat');     % load the packaged inputs

run('DualMarkerMainCode.m');          % produces variable "geneRANK"
head(geneRANK);                       % [gene, score] sorted descending
```

**2) 5-fold CV (50 repeats) at Top-100/150/200**

```matlab
% Requires "geneRANK" from the step above
run('KFoldCrossValidation.m');
```

**(Optional) Save the AUC lists**

```matlab
writematrix(GSE2034_aucScore_list_100, 'GSE2034_top100_auc.csv');
writematrix(GSE2034_aucScore_list_150, 'GSE2034_top150_auc.csv');
writematrix(GSE2034_aucScore_list_200,  'GSE2034_top200_auc.csv');
```

---

### Evaluating DualMarker

* **Top-K feature selection** is handled inside `KFoldCrossValidation.m` via `geneRANK(1:K)`.
* **Classifier** : by default we use a random forest (`classRF_*`). To use MATLAB’s native model:

```matlab
  % Example replacement using TreeBagger
  model = TreeBagger(500, train_X', categorical(train_Y));
  [~, score] = predict(model, test_X');
  score_pos = score(:,2);
  auc = AUC(test_Y, score_pos');
```

* **Standardization** : the script applies `zscore` + `mapstd`, consistent with the paper; feel free to swap in your preferred preprocessing.

---

### Quick test

```matlab
addpath(genpath(pwd));
load('Data/DualMarker_demo.mat');

run('DualMarkerMainCode.m');      % produces geneRANK
run('KFoldCrossValidation.m');    % prints AUC results
```

---
