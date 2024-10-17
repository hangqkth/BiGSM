# BiGSM: Bayeisan inference of GRN via sparse modelling

#### This the code for BiGSM, including algorithm implementation, experiments, results analysis and visualizaiton, and datasets.

## BiGSM_matlab
#### This folder includes the algorithm implemenation, experiments, datasets and data generation on Matlab, with some visualization of the results.

bigsm.m -- Implementation of the BiGSM algorithm.

rum_bigsm.m -- An example of how to run BiGSM using GeneSPIDER toolbox.

benchmark_test.m -- Excecute benchmark test using GeneSPIDER toolbox.

test_dream3/4/5 -- Run benchmark test with corresponding DREAM dataset.

genBenchmark.m -- Run benchmark test on dataser from GRNbenchmark.org.

test_ecoli.m -- Run benchmark test with biological e.coli. data.

After each benchmark test, save the result in .mat file. The results in .mat format will be analysis using python code.

## analysis_python
This folder includes the results analyzing and visualization. 

rank_grnbenchmark.py -- Visualize and analyse the results from GRNbenchmark.org.

plot_posterior -- Visualize the inferred posterior distributions from BiGSM and compared the point-estimates to the true GRN.

plot_ecoli -- Visualize the results on e.coli. biological data

ridgeline.py -- Analyse the density of inferred full GRN using ridgeline plots.



