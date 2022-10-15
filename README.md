# Environmental-oDNA-retention
Code for the paper 'Cellular and environmental dynamics influence species-specific extents of organelle gene retention'.

This code simulates a simple ODE model for nuclear and organelle gene expression in response to varying environments.

The code is presented as two versions: the publication version as a Jupyter notebook and an alternative (not identical) version in C and R. The C/R version has a wrapper bash script `run.sh` that compiles the simulation code, runs it, and then runs R to generate various summary plots. The payload is `param-scans.c` (scan through parameter space), `time-traces.c` (produce time series for specific parameterisations), and two `plot-X` R scripts.

The Python code uses `numpy`, `scipy`, `pandas`, `mpmath`, `matplotlib`, `seaborn`, and `bisect`. The R code uses `ggplot2` and `gridExtra`.

The `gene-counts` folder contains scripts to plot oDNA gene counts by clade from NCBI data. This requires the code from https://github.com/StochasticBiology/odna-loss . The `downloadorganelles` module of that code should be run. The `count-summary.sh` code here then recapitulates the `processorganelles` and `manuallabel` modules, and calls an R script which uses those outputs to plot the summary figures.
