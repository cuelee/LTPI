# LTPI 

`LTPI` (**L**iability **T**hreshold-based **P**henotype **I**mputation) is a statistical framework for predicting genetic liability of diseases by leveraging related phenotypes. It integrates EHR data using the GHK algorithm and maximum likelihood estimation, with automated trait-selection to improve risk prediction accuracy.

<br><br>

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Input Arguments](#input-arguments)
   - [Core Arguments](#core-arguments)
   - [Covariance Matrix Arguments](#covariance-matrix-arguments)
   - [Optional Arguments](#optional-arguments)
   - [Output Control](#output-control)
5. [Input File Formats](#input-file-formats)
   - [Phenotype Matrix](#1-phenotype-matrix---bin---con)
   - [Covariance Matrix](#2-covariance-matrix---gencov---envcov)
   - [Prevalence File](#3-prevalence-file---prevalence)
6. [Analysis Modes](#analysis-modes)
   - [Binary Phenotype Analysis](#1-binary-phenotype-analysis---bin)
   - [Continuous Phenotype Analysis](#2-continuous-phenotype-analysis---con)
   - [Trait Selection](#3-trait-selection---pick)
7. [Examples](#examples)
   - [Example 1: Binary Phenotype Analysis](#example-1-binary-phenotype-analysis)
   - [Example 2: Continuous Phenotype Analysis](#example-2-continuous-phenotype-analysis)
   - [Example 3: Trait Selection](#example-3-trait-selection)
8. [Contact](#contact)


<br><br>

## Overview
`LTPI` offers three analysis modes:
- **Binary phenotype analysis** (`--bin`)
- **Continuous phenotype analysis** (`--con`)
- **Trait selection optimization** (`--pick`)

### Key Features:
- Estimation of genetic liability using the GHK algorithm and maximum likelihood
- Covariance shrinkage to enhance model performance
- Automatic trait selection based on the `r2_o` criterion

<br><br>

## Installation
To install `LTPI`, clone the repository and ensure all dependencies are installed:
```bash
git clone https://github.com/your-repo/LTPI.git
cd LTPI
```

<br><br>

## Usage

`LTPI` can be run in three main modes depending on your input data:

1. **Binary phenotype analysis** (`--bin`)
2. **Continuous phenotype analysis** (`--con`)
3. **Trait selection analysis using ATSA** (`--pick`)

Provide the required input arguments based on your analysis type. Optional parameters, such as covariance matrices and analysis-specific parameters, can further adjust your analysis.

<br><br>

## Input Arguments

### Core Arguments

#### Binary Phenotype Analysis Input
| Argument       | Description                                             | Required?  |
|----------------|---------------------------------------------------------|------------|
| `--bin`        | Path to the binary phenotype input file.                 | Yes        |
| `--prevalence` | Path to disease prevalence information.                  | Yes        |
| `--gencov`     | Path to genetic covariance matrix.                       | Yes        |
| `--envcov`     | Path to environmental covariance matrix.                 | Optional   |

**Note**: The `--bin` can be used concurrently with `--pick` to perform trait selection optimization during phenotype analysis.

#### Continuous Phenotype Analysis Input
| Argument       | Description                                             | Required?  |
|----------------|---------------------------------------------------------|------------|
| `--con`        | Path to the quantitative phenotype input file.           | Yes        |
| `--bout`       | Same as the `--out` argument used in the `--bin` analysis. | Yes     |
| `--gencov`     | Path to genetic covariance matrix.                       | Yes        |
| `--envcov`     | Path to environmental covariance matrix.                 | Optional   |

**Note**: The `--con` can be used concurrently with `--pick` to perform trait selection optimization during phenotype analysis.

#### ATSA analysis Option (Optimization with `--pick`)
| Argument       | Description                                               | Required?  |
|----------------|-----------------------------------------------------------|------------|
| `--pick`       | Optimize non-target trait selection based on `r2_o`.      | Yes         |
| `--gencov`     | Path to genetic covariance matrix.                        | Yes        |
| `--pi`         | Target column name for trait selection (required for `--pick`). | Yes    |
| `--Q`          | Number of non-target traits to select (default: 30).      | No         |

### Covariance Matrix Arguments
| Argument             | Description                                                     | Required? |
|----------------------|-----------------------------------------------------------------|-----------|
| `--gencov`           | Path to genetic covariance matrix.                              | Yes       |
| `--envcov`           | Path to environmental covariance matrix.                        | Optional  |
| `--shrink`           | Apply covariance shrinkage: G for `gencov`, E for `envcov`, B for both. | Optional  |
| `--shrink_target`    | Target condition number for covariance shrinkage (default: 1.5).| No        |

**Note**: When `--envcov` is not provided, `LTPI` assumes the environmental covariance matrix (`envcov`) by transforming the genetic covariance matrix (`gencov`).

**Why Covariance Shrinkage?**
Covariance shrinkage is often necessary because covariance matrices can be non-positive semi-definite (non-PSD) due to bias, error, or subtle correlations at the boundary between PSD and non-PSD among pairs of multiple traits. This issue is more pronounced when working with a large number of traits, which reduces the stability of the LTPI analysis. Shrinking the covariance matrix ensures that the condition number is within a stable range, improving analysis reliability. A shrinkage with a condition number less than 10 is recommended for optimal performance.

<br><br>

### Optional Arguments
| Argument         | Description                                                    | Default  |
|-------------------------|----------------------------------------------------------------|------------|
| `--rint`         | Apply rank-based inverse normal transformation on LTPI scores. | False         |
| `--nsample`      | Number of samples for the GHK algorithm.                       | 50000      |
| `--r2`           | Threshold value for ATSA analysis.                             | 0.0        |

<br><br>

### Output Control
| Argument       | Description                                           | Default   |
|----------------|-------------------------------------------------------|-----------|
| `--out`        | Prefix for output files.                              | LTPI      |

<br><br>

## Input File Formats

### 1. Phenotype Matrix (`--bin`, `--con`)
A tab-delimited text file with the first column labeled `IID` for individual IDs. Subsequent columns contain trait names as headers and the corresponding phenotypic values for each individual.

### 2. Covariance Matrix (`--gencov`, `--envcov`)
A tab-delimited text file where the first row lists trait IDs. The remaining rows and columns contain the covariance values, forming an `n x n` matrix for `n` traits, with `n + 1` rows (including the header).

### 3. Prevalence File (`--prevalence`)
A two-column tab-delimited file. The first column, `TID`, contains trait IDs, and the second column, `prev`, contains prevalence values.

<br><br>

## Analysis Modes
`LTPI` provides three analysis modes based on input data type:

### 1. Binary Phenotype Analysis (`--bin`)
For binary phenotype data, this mode requires the binary input file (`--bin`), disease prevalence data (`--prevalence`), and covariance matrices (`--gencov`, `--envcov`). It estimates genetic liability for binary traits.

**Note**: At least one binary trait is required for analysis.

### 2. Continuous Phenotype Analysis (`--con`)
This mode handles continuous phenotype data. It requires the output from a previous binary analysis (`--bout`) and estimates liability based on continuous traits, integrating both binary and continuous data.

### 3. Trait Selection (`--pick`)
Optimizes non-target trait selection using the `r2_o` criterion, identifying traits most relevant to the target phenotype.

<br><br>

## Examples

The `LTPI` package includes a runnable example script, `runexample.bash`, which demonstrates how to use three analysis modes provided by `LTPI`. 

### Running the Example

To run the example, simply execute the `runexample.bash` script in the terminal:
```bash
bash runexample.bash
```

This will run all three modes of analysis (binary, continuous, and trait selection) using example data included in the package.

### Example 1: Binary Phenotype Analysis
```bash
python LTPI.py --bin .example/binary_input.txt --gencov ./example/genetic_covariance_bin.txt --prevalence ./example/prevalence.txt --out ./results/test_bin
```

### Example 2: Continuous Phenotype Analysis
```bash
python LTPI.py --con ./example/quantitative_input.txt --gencov ./example/genetic_covariance_con.txt --bout ./results/test_bin --out ./results/test_con
```

### Example 3: Trait Selection 
```bash
python LTPI.py --pick --pi trait_C --gencov ./example/genetic_covariance_bin.txt --out ./results/test_pick
```

<br><br>

## Contact
For support, inquiries, or further information, please submit an issue on the GitHub or contact the developer:

**Cue Hyunkyu Lee**  
Columbia University  
Email: [hl3565@cumc.columbia.edu](mailto:hl3565@cumc.columbia.edu?subject=[GitHub]%20LTPI%20Inquiry)