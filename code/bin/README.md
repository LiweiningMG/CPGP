# Executable programs

This directory stores executable programs required to reproduce the analyses in this repository. The workflow may call the following programs:

```text
plink          PLINK, used for genotype processing, PCA and LD calculations
mtg2           MTG2, used to construct multi-population genomic relationship matrices
dmuai          DMUAI, used for variance component estimation
dmu1           DMU1, used for GBLUP prediction
run_dmuai      Helper script for running DMUAI
run_dmu4       Helper script for running DMU
gmatrix        Helper program for genomic relationship matrix processing
ldblock        Helper program for LD block definition
LD_mean_r2     Helper program for summarising LD decay
mbBayesABLD    Multi-population Bayesian genomic prediction program used in this study
QMSim_selected Auxiliary executable retained for compatibility with the original workflow
```

These programs are provided only to reproduce the analyses in this project. Any other use must comply with the licence terms, user agreements or distribution policies issued by the original software authors, development teams or institutions. If equivalent programs are already installed in the system environment, users may edit the executable paths in `main.sh` or related shell scripts to call the system versions instead.

Before running the analysis, make sure the programs in this directory have executable permissions:

```bash
chmod +x code/bin/*
```

The Chinese version of this note is available in `readme_cn.md`.
