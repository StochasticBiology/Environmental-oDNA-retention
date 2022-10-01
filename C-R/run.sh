# compile and run code to generate specific time series for some experiments (should be fast)
gcc -o3 time-traces.c -lm -o time-traces.ce
./time-traces.ce > time-traces-out.tmp

# compile and run code to scan through parameter spaces for various experiments (will take longer)
gcc -o3 param-scans.c -lm -o param-scans.ce
./param-scans.ce > param-scans-out.tmp

# produce plots of time trace experiments
Rscript plot-time-traces.R

# produce plots of parameter scan experiments
Rscript plot-param-scans.R
