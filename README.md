# LTPI 

`LTPI` is a command line tool to derive new phenotypes for a target disease by combining genetic relatedness information with phenotypic information.

## Getting Started

### Installation

To install `LTPI`, clone this repository using the following commands:

```bash
git clone https://github.com/cuelee/LTPI.git
cd LTPI
```

### Dependencies

`LTPI` requires the following Python dependencies to be installed:

- Python 3.x
- pandas
- numpy
- numba
- scipy
- scikit-learn

You can install these dependencies using pip or conda.

### Creating the LTPI Environment with Conda

If users prefer to use conda, they can define the LTPI environment using the following command:

```bash
conda env create -f environment.yml
conda activate LTPI
```

## Running LTPI

After installing the dependencies, you can start using `LTPI` by running:

```bash
./LTPI.py -h
```

This command displays the help information and available commands.

## Updating LTPI

To update `LTPI` to the latest version, navigate to your `LTPI/` directory and run:

```bash
git pull
```

If `LTPI` is already up to date, the output will be:

```
Already up to date.
```

## For First-Time Users

New users can find an executable example file in the `LTPI/` directory, named `runexample.bash`. This example provides a practical guide to using `LTPI`.

## Citation

If you utilize `LTPI` in your research, please cite the following paper:

LEE, C. H., Khan, A., Weng, C., Buxbaum, J. D., Kiryluk, K., Ionita-Laza, I., "Liability threshold model-based disease risk prediction based on Electronic Health Record phenotypes," Under Review.

## Support

For any issues or queries regarding `LTPI`, feel free to reach out via email: hl3565@cumc.columbia.edu

## License 

`LTPI` is released under the MIT License. See the LICENSE file for more details.

## Authors

- Cue Hyunkyu Lee, Columbia University Irving Medical Center

---

