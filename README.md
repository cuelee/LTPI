# LTPI

`LTPI` (**L**iability **T**hreshold-based **P**henotype **I**ntegration) is a statistical framework designed to predict genetic liability for diseases by leveraging related phenotypes. It uses EHR data with the GHK algorithm and maximum likelihood estimation, featuring automated trait-selection to enhance risk prediction accuracy.

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
   - [Prevalence File](#prevalence-file)
6. [Analysis Modes](#analysis-modes)
7. [Examples](#examples)
8. [Contact](#contact)

---

## Overview
`LTPI` operates in three main analysis modes:
- **Binary phenotype analysis** (`--bin`)
- **Continuous phenotype analysis** (`--con`)
- **Trait selection optimization** (`--pick`)

### Key Features:
- Predicts genetic liability using the GHK algorithm and maximum likelihood estimation.
- Incorporates covariance shrinkage to improve model stability.
- Implements automatic trait selection using the `r2_o` criterion.

---

## Installation

### Option 1: Install via `environment.yml`

To install `LTPI` using a conda environment, follow these steps:

1. Clone the repository:
    ```bash
    git clone https://github.com/your-repo/LTPI.git
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
    git clone https://github.com/your-repo/LTPI.git
    cd LTPI
    ```

2. Install the required dependencies manually using `pip`:
    ```bash
    pip install numpy pandas scipy scikit-learn numba argparse
    ```
---

## Usage

`LTPI` offers three analysis modes based on the input data:

1. **Binary Phenotype Analysis (`--bin`)**: Estimates genetic liability for binary traits. Requires binary input, prevalence data, and covariance matrices.
2. **Continuous Phenotype Analysis (`--con`)**: Handles continuous traits and integrates them with binary traits. Requires previous binary output (`--bout`).
3. **Trait Selection Optimization (`--pick`)**: Automatically selects non-target traits to optimize genetic liability prediction.

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
| `--pick`       | Optimize non-target trait selection based on `r2_o`.      | Yes         |
| `--gencov`     | Path to genetic covariance matrix.                        | Yes        |
| `--pi`         | Target column name for trait selection.                   | Yes        |
| `--Q`          | Number of non-target traits to select (default: 30).      | No         |

### Covariance Matrix Arguments
| Argument          | Description                                                     | Required? |
|-------------------|-----------------------------------------------------------------|-----------|
| `--gencov`        | Path to genetic covariance matrix.                              | Yes       |
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

### Phenotype Matrix (`--bin`, `--con`)
A tab-delimited text file. The first column must be labeled `IID` for individual IDs, followed by columns for each trait.

### Covariance Matrix (`--gencov`, `--envcov`)
A tab-delimited matrix file with the first row containing trait IDs and subsequent rows forming an `n x n` matrix of covariance values for `n` traits.

### Prevalence File (`--prevalence`)
A two-column tab-delimited file with `TID` for trait IDs and `prev` for disease prevalence. This file is required for binary traits.

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
