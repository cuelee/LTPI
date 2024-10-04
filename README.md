# LTPI: Liability Threshold Model for Phenotype Imputation

`LTPI` is a statistical framework for predicting the genetic liability of target diseases by leveraging related phenotypes. It integrates phenotypic data from electronic health records (EHR) using the GHK algorithm and maximum likelihood estimation. LTPI’s automated trait-selection algorithm enhances the accuracy of disease risk prediction.
<br><br>
## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Input Arguments](#input-arguments)
5. [Analysis Modes](#analysis-modes)
6. [Examples](#examples)
<br><br>
## Overview
`LTPI` provides three different analysis modes:
- **Binary phenotype analysis** (`--bin`)
- **Continuous phenotype analysis** (`--con`)
- **Trait selection optimization** (`--pick`)

### Included Features:
- Evaluation of expected genetic liability using GHK algorithm and maximum likelihood
- Covariance shrinkage to improve the model's performance
- Automatic trait selection based on `r2_o` criteria
<br><br>
## Installation
To install `LTPI`, clone the repository and ensure all dependencies are installed:
```bash
git clone https://github.com/your-repo/LTPI.git
cd LTPI
```
<br><br>
## Usage
Run `LTPI` with one of the available analysis modes by providing the required input arguments for your analysis. There are three main analysis modes:
1. **Binary phenotype input (`--bin`)**
2. **Continuous phenotype input (`--con`)**
3. **Trait selection optimization (`--pick`)**

You can customize the analysis by using various optional parameters, such as specifying covariance matrices, optimizing trait selection, and applying covariance shrinkage.
<br><br>
## Input Arguments
### Output Arguments
| Argument       | Description                                | Default   |
|----------------|--------------------------------------------|-----------|
| `--out`        | Output file prefix.                        | `LTPI`    |

### Test Mode Arguments
#### Binary Phenotype Input
| Argument       | Description                                             | Required?  |
|----------------|---------------------------------------------------------|------------|
| `--bin`        | Path to the binary phenotype input file.                 | Yes        |
| `--prevalence` | Path to disease prevalence information.                  | Yes        |
| `--gencov`     | Path to genetic covariance matrix.                       | Yes        |
| `--envcov`     | Path to environmental covariance matrix.                 | Optional   |

#### Continuous Phenotype Input
| Argument       | Description                                             | Required?  |
|----------------|---------------------------------------------------------|------------|
| `--con`        | Path to the quantitative phenotype input file.           | Yes        |
| `--bout`       | Same as the `--out` argument used in the `--bin` test.   | Yes        |
| `--gencov`     | Path to genetic covariance matrix.                       | Yes        |
| `--envcov`     | Path to environmental covariance matrix.                 | Optional   |

#### Trait Selection (Optimization with `--pick`)
| Argument       | Description                                               | Required?  |
|----------------|-----------------------------------------------------------|------------|
| `--pick`       | Optimize non-target trait selection based on `r2_o`.       | No         |
| `--pi`         | Target column name for trait selection (required for `--pick`). | Yes    |
| `--Q`          | Number of non-target traits to select (default: 30).       | No         |

### Rank-Based Inverse Normal Transformation
| Argument       | Description                                               | Required?  |
|----------------|-----------------------------------------------------------|------------|
| `--rint`       | Apply rank-based inverse normal transformation on LTPI scores. | No     |

### Covariance Matrix Arguments
| Argument          | Description                                                     | Required? |
|-------------------|-----------------------------------------------------------------|-----------|
| `--gencov`        | Path to genetic covariance matrix.                              | Yes       |
| `--envcov`        | Path to environmental covariance matrix.                        | Optional  |
| `--shrink`        | Apply covariance shrinkage: G for GENCOV, E for ENVCOV, B for both. | Optional  |
| `--shrink_target` | Target condition number for covariance shrinkage (default: 1.5).| No        |

**Note**: When `--envcov` is not provided, `LTPI` assumes the environmental covariance matrix (`envcov`) by transforming the genetic covariance matrix (`gencov`).

### Parameters Specific to GHK Algorithm
| Argument           | Description                                       | Default |
|--------------------|---------------------------------------------------|---------|
| `--nsample_main`   | Number of samples for the GHK algorithm.          | 50000   |

### Parameters Specific to R2 Selection (ATSA)
| Argument   | Description                                    | Default |
|------------|------------------------------------------------|---------|
| `--r2`     | `r2_o` threshold for ATSA.                     | 0.0     |
<br><br>
## Analysis Modes
`LTPI` offers three distinct modes of analysis, depending on the input data type:

### 1. Binary Phenotype Input (`--bin`)
This mode is designed for binary phenotype data. It requires the binary phenotype input file using `--bin`, along with disease prevalence information (`--prevalence`) and covariance matrices (`--gencov` and `--envcov`). The results will estimate the genetic liability for binary traits.

**Note**: At least one binary target trait must be provided. The `LTPI` analysis cannot proceed if the number of binary traits is zero.

### 2. Continuous Phenotype Input (`--con`)
This mode handles continuous phenotype data. It requires a reference to a previous binary phenotype analysis (using --bout) and performs an additional estimation based on continuous traits, leveraging both the binary and continuous information.

### 3. Trait Selection (`--pick`)
In this mode, `LTPI` optimizes the selection of non-target traits based on the `r2_o` criterion. It is useful for identifying key features that are most relevant to the target phenotype.
<br><br>
## Examples

The `LTPI` package includes a runnable example script called `runexample.bash`, which demonstrates how to use the different analysis modes provided by `LTPI`. This script allows you to easily test the package and explore its functionality.

### Running the Example

To run the example, simply execute the `runexample.bash` script in the terminal:
```bash
bash runexample.bash
```

This will run all three modes of analysis (binary, continuous, and trait selection) using example data included in the package.

### Example 1: Binary Phenotype Analysis
```bash
python LTPI.py --bin ./example/binary_input.txt --gencov ./example/genetic_covariance_bin.txt --prevalence ./example/prevalence.txt --out ./results/test_bin
```

### Example 2: Continuous Phenotype Analysis
```bash
python LTPI.py --con ./example/quantitative_input.txt --gencov ./example/genetic_covariance_con.txt --bout ./results/test_bin --out ./results/test_con
```

### Example 3: Trait Selection Optimization
```bash
python LTPI.py --pick --pi trait_C --gencov ./example/genetic_covariance_con.txt --out ./results/test_pick
```
<br><br>
## Contact
For any issues or questions, feel free to raise an issue on GitHub or contact the author.
