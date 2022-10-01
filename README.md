# Environmental-oDNA-retention
Code for the paper 'Cellular and environmental dynamics influence species-specific extents of organelle gene retention'.

This code simulates a simple ODE model for nuclear and organelle gene expression in response to varying environments.

The code is presented as two versions: the publication version as a Jupyter notebook and an alternative (not identical) version in C and R. The C/R version has a wrapper bash script `run.sh` that compiles the simulation code, runs it, and then runs R to generate various summary plots. 

The Python code uses `numpy`, `scipy`, `pandas`, `mpmath`, `matplotlib`, `seaborn`, and `bisect`. The R code uses `ggplot2` and `gridExtra`.
