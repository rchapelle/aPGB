# Parallel gradient boosting for joint density estimation

This repository provides the source code of the experiments for the paper: [Nonparametric autoregressive estimation of probability distributions with parallel gradient boosting](https://arxiv.org/abs/2610.31524).

## Install the R package

Our method is available within the R package [opencesp](https://cran.r-project.org/web/packages/opencesp/index.html) available on CRAN. Easily install it with:  

```
install.packages("opencesp")
```

## Experiments

To run the experiments, first clone the repository with:
```
git clone https://github.com/rchapelle/aPGB.git
```
Then, install the required R packages listed at the beginning of the `script.R` file, e.g.:
```
install.packages("np")
```
Finally, in the `script.R` file, set the variable `wd` to the directory in which the repository was cloned, and `plot_dir` to the directory in which you want the plots to be saved. All experiments are contained in the `script.R` file and can then be conducted by running the script, for example within RStudio.

## Datasets

The [Hepatitis dataset](https://doi.org/10.24432/C5D612) is licensed under the Creative Commons Attribution 4.0 International license.

## Reference

If you use this repository in your research, please cite our work using the following reference:

```
@misc{chapelle2026nonparametricautoregressiveestimationprobability,
      title={Nonparametric autoregressive estimation of probability distributions with parallel gradient boosting}, 
      author={Rémy Chapelle},
      year={2026},
      eprint={2610.31524},
      archivePrefix={arXiv},
      primaryClass={stat.ML},
      url={https://arxiv.org/abs/2610.31524}, 
}
```
