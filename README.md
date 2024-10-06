# LTPI

`LTPI` (**L**iability **T**hreshold-based **P**henotype **I**ntegration) is a statistical framework designed to predict genetic liability for diseases by leveraging related phenotypes in EHR. The risk score construction utilizes the GHK algorithm and maximum likelihood estimation, with efficacy that can be further enhanced through automated feature selection.

---

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
   - [Phenotype Matrix](#phenotype-matrix)
   - [Covariance Matrix](#covariance-matrix)
   - [Prevalence](#prevalence-of-binary-traits)
6. [Output Files](#output-files)
7. [Examples](#examples)
8. [Contact](#contact)

---

## Overview

`LTPI` operates through three analysis modes:

1. **Binary phenotype analysis** (`--bin`): Merges all binary traits to estimate genetic liability.
2. **Continuous phenotype analysis** (`--con`): Combines the results from `--bin` with continuous traits to generate the final LTPI score.
3. **Trait selection optimization** (`--pick`): Identifies and selects the most important traits based on the `r2_o` criterion, refining the analysis for optimal results.

**Important:** Two steps—`--bin` and `--con`—must be run sequentially to obtain the final LTPI score. Skipping any step may result in incomplete analysis. `--pick` is optional, but it usually improves the performance.

### Key Features:

- **GHK Algorithm and Maximum Likelihood Estimation**: Provides accurate genetic liability estimates.
- **Covariance Shrinkage**: Enhances model stability by applying shrinkage to covariance matrices, especially when analyzing a large number of traits.
- **Automated Trait Selection**: Utilizes the `r2_o` criterion to automatically select the most relevant non-target traits, enhancing prediction accuracy.

---

## Installation

### Option 1: Install via `environment.yml`

To install `LTPI` using a conda environment, follow these steps:

1. Clone the repository:
    ```bash
    git clone https://github.com/cuelee/LTPI.git
    cd LTPI
    ```

2. Create the conda environment from the `environment.yml` file:
    ```bash
    conda env create -f environment.yml
    ```

3. Activate the environment:
    ```bash
    conda activate LTPI
    ```

### Option 2: Manual Installation

If you prefer to manually install the dependencies, follow these steps:

1. Clone the repository:
    ```bash
    git clone https://github.com/cuelee/LTPI.git
    cd LTPI
    ```

2. Install the required dependencies manually using `pip`:
    ```bash
    pip install numpy pandas scipy scikit-learn numba argparse
    ```
---

## Usage

`LTPI` provides three analysis modes depending on the input data:

1. **Binary Phenotype Analysis (`--bin`)**: Estimates genetic liability based on binary traits. Requires binary input, prevalence data, and covariance matrices (`--gencov`, `--envcov`).
2. **Continuous Phenotype Analysis (`--con`)**: Integrates continuous traits with the results from binary trait analysis. Requires the binary output from `--bin` (`--bout`), continuous input, and covariance matrices.
3. **Trait Selection Optimization (`--pick`)**: Optionally selects the most relevant non-target traits to enhance genetic liability prediction. Requires a disease of interest (`--pi`) and genetic covariance matrices (`--gencov`), applicable to either `--bin` or `--con`.

Specify the required core arguments for each mode and optional parameters like covariance matrices to fine-tune the analysis.

---

## Input Arguments

### Core Arguments

#### Binary Phenotype Analysis
| Argument       | Description                                             | Required?  |
|----------------|---------------------------------------------------------|------------|
| `--bin`        | Path to binary phenotype input file.                     | Yes        |
| `--prevalence` | Path to disease prevalence file.                         | Yes        |
| `--gencov`     | Path to genetic covariance matrix.                       | Yes        |
| `--envcov`     | Path to environmental covariance matrix.                 | Optional   |

#### Continuous Phenotype Analysis
| Argument       | Description                                             | Required?  |
|----------------|---------------------------------------------------------|------------|
| `--con`        | Path to continuous phenotype input file.                 | Yes        |
| `--bout`       | Output file from previous binary analysis.               | Yes        |
| `--gencov`     | Path to genetic covariance matrix.                       | Yes        |
| `--envcov`     | Path to environmental covariance matrix.                 | Optional   |

#### Trait Selection (Optimization with `--pick`)
| Argument       | Description                                               | Required?  |
|----------------|-----------------------------------------------------------|------------|
| `--pick`       | Optimize non-target trait selection based on `r2_o`.      | Yes        |
| `--gencov`     | Path to genetic covariance matrix.                        | Yes        |
| `--pi`         | Target column name for trait selection.                   | Yes        |
| `--Q`          | Number of non-target traits to select (default: 30).      | No         |

### Covariance Matrix Arguments
| Argument          | Description                                                     | Required? |
|-------------------|-----------------------------------------------------------------|-----------|
| `--gencov`        | Path to genetic covariance matrix.                              | Optional       |
| `--envcov`        | Path to environmental covariance matrix.                        | Optional  |
| `--shrink`        | Apply covariance shrinkage: G for `gencov`, E for `envcov`, B for both. | Optional  |
| `--shrink_target` | Target condition number for covariance shrinkage (default: 1.5). | No        |

### Optional Arguments
| Argument    | Description                                                    | Default |
|-------------|----------------------------------------------------------------|---------|
| `--rint`    | Apply rank-based inverse normal transformation to LTPI scores. | False   |
| `--nsample` | Number of samples for GHK algorithm.                           | 50000   |
| `--r2`      | Threshold value for trait selection analysis.                  | 0.0     |

### Output Control
| Argument | Description                  | Default |
|----------|------------------------------|---------|
| `--out`  | Prefix for output files.     | LTPI    |

---

## Input File Formats

### Phenotype Matrix
(Used with --bin, --con)
A tab-delimited text file. The first column must be labeled `IID` for individual IDs, followed by columns for each trait.

### Covariance Matrix
(Used with`--gencov`, `--envcov`)
A tab-delimited matrix file with the first row containing trait IDs and subsequent rows forming an `n x n` matrix of covariance values for `n` traits.

### Prevalence of binary traits
(Used with `--prevalence`)
A two-column tab-delimited file with `TID` for trait IDs and `prev` for disease prevalence. This file is required for binary traits.

---

## Output Files

After running the LTPI analysis, the following output files will be generated with the corresponding suffixes:

1. **{out}.config_summary**  
   Contains three columns:
   - `CONF`: Configurations (case-control status).
   - `EG`: Estimated genetic liabilities.
   - `G_SE`: Standard error of EG.

2. **{out}.ltpi_scores.gz**  
   Contains three columns:
   - `IID`: Individual ID.
   - `LTPI`: LTPI score for each individual.
   - `CONF`: Configurations (case-control status).

3. **{out}.shrinkage_value**  
   Contains a single value representing the shrinkage value used in the ATSA analysis.

4. **{out}.r2_summary**  
   A three-column file with multiple rows, where each row represents a trait:
   - `TID`: Trait ID
   - `R2`: Cumulative R2 for the trait and all traits listed above it.
   - `COND`: Estimated condition number for the genetic covariance matrix corresponding to the cumulative R2.

6. **{out}.selected_traits**  
   A single-column file listing the traits selected by ATSA.

---

## Examples

`LTPI` includes a runnable example script (`runexample.bash`). To run all analysis modes:

```bash
bash runexample.bash
```

### Example 1: Binary Phenotype Analysis
```bash
python LTPI.py --bin example/binput.txt --gencov example/genetic_covariance_bin.txt --prevalence example/prevalence.txt --out test_bin
```

### Example 2: Continuous Phenotype Analysis
```bash
python LTPI.py --con example/qinput.txt --gencov example/genetic_covariance_con.txt --bout test_bin --out test_con
```

### Example 3: Trait Selection
```bash
python LTPI.py --pick --pi trait_C --gencov example/genetic_covariance_bin.txt --out test_pick
```

---

## Contact

For inquiries, support, or feedback, submit an issue on GitHub or contact:

**Cue Hyunkyu Lee**  
Columbia University  
Email: [hl3565@cumc.columbia.edu](mailto:hl3565@cumc.columbia.edu?subject=[GitHub]%20LTPI%20Inquiry)
