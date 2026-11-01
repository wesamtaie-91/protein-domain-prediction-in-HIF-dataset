# protein-domain-prediction-in-HIF-dataset

# Protein Domain Prediction and Protein–Protein Interaction Analysis Using CDM

## Overview

This repository contains a MATLAB implementation of a computational framework for **protein domain prediction and protein–protein interaction (PPI) analysis** using the Common Domain Model (CDM).

The proposed framework identifies informative sequence features using overlapping 3-mers and Laplacian Score, groups proteins using K-means clustering, infers missing domain annotations, and predicts potential protein–protein interactions based on Jaccard similarity between protein domain sets.

The implementation is organized into two MATLAB scripts corresponding to the main stages of the proposed algorithm.

## Repository Structure

```text
.
├── Protein_Feature_Clustering_Steps1to6.m
├── Protein_Domain_CDM_Steps7to10.m
├── ProteinData.xlsx
├── PPIData.xlsx
└── README.md
```

## Requirements

* MATLAB R2022a or later.
* Statistics and Machine Learning Toolbox.
* Input data in Microsoft Excel (`.xlsx`) format.

## Algorithm Implementation

### 1. Feature Extraction and Clustering

**Script:** `Protein_Feature_Clustering_Steps1to6.m`

This script implements Steps 1–6 of the proposed framework:

1. **Data Integration:** Integrates protein sequences with available Pfam domain annotations.
2. **Data Cleaning:** Removes duplicate sequences, redundant records, and invalid sequences.
3. **Feature Extraction:** Generates overlapping 3-mer frequency vectors from protein sequences.
4. **Feature Scoring:** Calculates Laplacian Scores to assess the relevance of extracted features.
5. **Feature Selection:** Selects the most informative features.
6. **K-means Clustering:** Groups proteins into eight clusters based on their selected feature vectors.

**Output files:**

* `ProteinClusteringResults.mat`: Stores processed protein data, selected features, Laplacian Scores, and cluster assignments.
* `SelectedFeatures.xlsx`: Contains the selected sequence features.

### 2. Domain Completion and PPI Prediction

**Script:** `Protein_Domain_CDM_Steps7to10.m`

This script implements Steps 7–10 of the proposed framework:

7. **Domain Completion:** Infers missing domain annotations using the dominant known domains within each cluster.
8. **Domain Prediction Evaluation:** Compares predicted domains with available known annotations.
9. **PPI Prediction:** Estimates protein–protein interactions using Jaccard similarity between protein domain sets.
10. **PPI Evaluation:** Calculates performance metrics using available reference interaction labels.

**Output files:**

* `PredictedProteinDomains.xlsx`: Contains predicted protein domain annotations and confidence scores.
* `PredictedPPI.xlsx`: Contains PPI prediction scores and classifications.
* `PPI_Evaluation_Metrics.xlsx`: Contains the calculated PPI evaluation metrics.

## Input Data

The implementation requires two Excel input files.

### ProteinData.xlsx

Contains protein sequences and their available Pfam annotations.

| Column        | Description                   |
| ------------- | ----------------------------- |
| `ProteinID`   | Unique protein identifier     |
| `Sequence`    | Amino acid sequence           |
| `PfamDomains` | Known Pfam domain identifiers |

Multiple domain identifiers should be separated by semicolons (`;`). An empty field indicates that no domain annotation is available.

**Example:**

| ProteinID | Sequence | PfamDomains     |
| --------- | -------- | --------------- |
| P001      | MKK...   | PF00001;PF00002 |
| P002      | MAA...   | PF00003         |
| P003      | MRL...   |                 |

### PPIData.xlsx

Contains reference protein–protein interaction pairs.

| Column        | Description                                                       |
| ------------- | ----------------------------------------------------------------- |
| `ProteinA`    | Identifier of the first protein                                   |
| `ProteinB`    | Identifier of the second protein                                  |
| `Interaction` | Reference interaction label (1 = interaction, 0 = no interaction) |

**Example:**

| ProteinA | ProteinB | Interaction |
| -------- | -------- | ----------- |
| P001     | P002     | 1           |
| P001     | P003     | 0           |

**Note:** Protein identifiers must match those in `ProteinData.xlsx`.

## Algorithm Parameters

The current implementation uses the following parameters:

| Parameter                    | Value |
| ---------------------------- | ----: |
| k-mer length                 |     3 |
| Amino acid alphabet size     |    20 |
| Possible 3-mer features      | 8,000 |
| Selected features            |    92 |
| Number of clusters           |     8 |
| K-means replicates           |     5 |
| Maximum iterations           |   300 |
| Nearest neighbors            |    10 |
| Jaccard similarity threshold |   0.3 |

These parameters can be modified directly in the MATLAB scripts.

## How to Run

1. Clone or download this repository.

2. Place `ProteinData.xlsx` and `PPIData.xlsx` in the same directory as the MATLAB scripts.

3. Open MATLAB.

4. Set the current working directory to the repository folder.

5. Run the first script:

   ```matlab
   Protein_Feature_Clustering_Steps1to6
   ```

6. After successful execution, the script generates `ProteinClusteringResults.mat`.

7. Run the second script:

   ```matlab
   Protein_Domain_CDM_Steps7to10
   ```

8. The resulting domain predictions, PPI predictions, and evaluation metrics will be saved as Excel files in the working directory.

**Important:** The second script depends on the output generated by the first script. Therefore, both scripts must be executed in the specified order.

## Methodological Notes

* The implementation uses the standard Laplacian Score, where lower scores indicate more relevant features. Accordingly, the current script selects the 92 features with the lowest scores.
* The 3-mer representation considers the 20 standard amino acids, resulting in 8,000 possible 3-mer features.
* Domain completion currently uses the dominant known domain within each cluster. This is a simplified implementation and does not fully account for proteins containing multiple domains.
* The current domain evaluation is descriptive. An independent test set with hidden domain annotations is required for unbiased evaluation.
* PPI prediction is based on Jaccard similarity between protein domain sets. Additional domain–domain interaction rules may be incorporated depending on the intended CDM formulation.
* PPI evaluation considers only protein pairs with available reference labels and valid prediction scores.

## Computational Considerations

The proposed framework is designed for protein sequence datasets. However, processing hundreds of thousands of sequences can require substantial computational resources, particularly during similarity-graph construction and nearest-neighbor searches.

For large-scale datasets, memory optimization, batch processing, or approximate nearest-neighbor techniques may be required.

## Workflow

```text
Protein Sequences and Pfam Annotations
                 |
                 v
          Data Integration
                 |
                 v
           Data Cleaning
                 |
                 v
       Overlapping 3-mer Extraction
                 |
                 v
        Laplacian Score Calculation
                 |
                 v
            Feature Selection
                 |
                 v
           K-means Clustering
                 |
                 v
           Domain Completion
                 |
                 v
       Domain Prediction Evaluation
                 |
                 v
          CDM-based PPI Prediction
                 |
                 v
           PPI Performance Evaluation
```

## License

The source code is made available for academic research and peer-review purposes.
All rights are reserved by the authors.

A formal open-source license will be specified upon publication of the associated research article.

## Citation

## Citation

If you use this implementation in your research, please cite the associated research article.
manuscript title "A Novel Approach for Protein Domain Prediction in Human Intestinal Flora:Integrating K-mer Analysis with K-means Clustering"

**Manuscript status:** Under review.

The complete citation and DOI will be added following publication.
