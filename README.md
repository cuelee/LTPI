# LTPI

`LTPI` (Liability Threshold-based Phenotype Integration) estimates the conditional expected genetic liability of a target disease using multi-trait phenotypic data and external covariance structure. The risk score construction utilizes the GHK algorithm and maximum likelihood estimation, with efficacy that can be further enhanced through automated feature selection.

---

> **Update (March 22, 2026)**  
> `--pi` is now a required argument that explicitly defines the target phenotype. It is no longer inferred from column order.  

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
   - [Trait Alignment](#trait-alignment)
5. [Input File Formats](#input-file-formats)
   - [Phenotype Matrix](#phenotype-matrix)
   - [Covariance Matrix](#covariance-matrix)
   - [Prevalence of binary traits](#prevalence-of-binary-traits)
6. [Output Files](#output-files)
7. [Examples](#examples)
8. [Citation](#citation)

---

## Overview

`LTPI` (Liability Threshold-based Phenotype Integration) estimates the conditional expected genetic liability of a target disease by integrating multiple phenotypes under a liability threshold model.

The method combines:
- binary and continuous phenotypes
- external genetic and environmental covariance structures
- GHK-based sampling and maximum likelihood estimation

`LTPI` operates in three modes:

1. **Binary analysis (`--bin`)**  
   Estimates liability using binary traits.

2. **Continuous analysis (`--con`)**  
   Extends binary results by incorporating continuous traits.

3. **Trait selection (`--pick`)**  
   Optionally selects informative traits using an R2-based criterion.

**Pipeline requirement:**
- `--bin` must be run before `--con`
- `--con` requires `--bout` from the binary step

**Important:**  
The target phenotype (`--pi`) is automatically placed as the first column across all inputs. All computations assume this ordering.

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
2. **Continuous Phenotype Analysis (`--con`)**:  Extends the binary analysis by incorporating continuous traits. The target phenotype (`--pi`) is obtained from the binary output (`--bout`) and is not required in the continuous phenotype input.
3. **Trait Selection Optimization (`--pick`)**: Optionally selects the most relevant non-target traits to enhance genetic liability prediction. Requires a disease of interest (`--pi`) and genetic covariance matrices (`--gencov`), applicable to either `--bin` or `--con`.

Specify the required core arguments for each mode and optional parameters like covariance matrices to fine-tune the analysis.

---

## Input Arguments

#### Core Arguments

| Argument       | Description | Required |
|----------------|-------------|----------|
| `--bin`        | Path to binary phenotype input file. | Yes (binary mode) |
| `--con`        | Path to continuous phenotype input file. | Yes (continuous mode) |
| `--prevalence` | Path to prevalence file for binary traits. | Yes (binary mode) |
| `--bout`       | Output prefix from the binary step. | Yes (continuous mode) |
| `--pi`         | Target phenotype used for prediction. | Yes |

### Covariance Matrix Arguments

| Argument          | Description | Required |
|-------------------|-------------|----------|
| `--gencov`        | Genetic covariance matrix. | Yes |
| `--envcov`        | Environmental covariance matrix (optional; defaults to 1 - diag(GENCOV)). | No |
| `--shrink`        | Covariance shrinkage mode: G, E, or B. | No |
| `--shrink_target` | Target condition number for shrinkage. | No |

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

### Trait Alignment

LTPI automatically aligns traits across phenotype inputs and covariance matrices.

- In binary mode (`--bin`), the target phenotype (`--pi`) must be present in the phenotype input.
- In continuous mode (`--con`), the target phenotype is obtained from the binary output (`--bout`) and does not need to appear in the continuous phenotype input.
- Only traits shared between inputs and covariance matrices are used in the analysis.

Incorrect or inconsistent trait naming across files may result in errors.

---

## Input File Formats

### Phenotype Matrix
(Used with `--bin`, `--con`)  
A tab-delimited text file where the first column must be labeled `IID` for individual IDs, followed by columns for each trait.

- **Binary traits** are coded as `1` (case) and `0` (control), with `NA` for missing values.
- **Continuous traits** should follow a standard normal distribution with a mean of `0` and variance of `1`. Missing phenotypes should be labeled as `NA`. If your data does not meet these criteria, you can use the **Rank-Based Inverse Normal Transformation (RINT)** for normalization. This is optional but recommended for continuous data.

You can find the `RINT.py` script for performing this transformation in the utility folder of the repository [LINK](https://github.com/cuelee/LTPI/blob/main/utility/RINT.py).


### Covariance Matrix
(Used with `--gencov`, `--envcov`)  
A tab-delimited matrix file where the first row contains trait IDs, and the subsequent rows form an `n x n` matrix of covariance values for `n` traits. Missing values (`NA`) are not allowed.

### Prevalence of binary traits
(Used with `--prevalence`)
A two-column tab-delimited file with `TID` for trait IDs and `prev` for disease prevalence. This file is required for binary traits.

**Note**: We assume that the biobank EHR data approximates population prevalence. However, users should consider potential biases introduced by the biobank's data collection process, especially if the data contains sampling bias.

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
python LTPI.py --bin example/binput.txt --pi trait_A --gencov example/genetic_covariance_bin.txt --prevalence example/prevalence.txt --out test_bin
```

### Example 2: Continuous Phenotype Analysis
```bash
python LTPI.py --con example/qinput.txt --pi trait_A --gencov example/genetic_covariance_con.txt --bout test_bin --out test_con
```

### Example 3: Trait Selection
```bash
python LTPI.py --pick --pi trait_C --gencov example/genetic_covariance_bin.txt --out test_pick
```

---

## Acknowledgments 
We thank Matthieu Pluntz for helpful comments and suggestions that improved the code.

---

## Contact

For inquiries, support, or feedback, submit an issue on GitHub or contact:

**Cue Hyunkyu Lee**  
Columbia University  
Email: [hl3565@cumc.columbia.edu](mailto:hl3565@cumc.columbia.edu?subject=[GitHub]%20LTPI%20Inquiry)

## Citation

If you use LTPI in your work, please cite:

Lee, C.H., Khan, A., Wang, C. et al.  
*Liability threshold model-based disease risk prediction based on electronic health record phenotypes.*  
Nature Genetics 57, 2872–2881 (2025).  
https://doi.org/10.1038/s41588-025-02370-4