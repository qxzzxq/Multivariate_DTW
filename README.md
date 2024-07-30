# Multivariate Time Series Clustering for HPC Performance Monitoring
This repository contains the MATLAB code for the paper Multivariate Time Series Clustering for HPC Performance Monitoring.

To run the code in MATLAB for the first time, follow these steps:

## Build the `dtwCPP` MEX Function:

In the MATLAB command window, execute the following command to compile the `dtwCPP` MEX function

```
mex dtwCPP.cpp
```

## Data Organization

Each folder should contain data files in `.mat` format, corresponding to the name of the folder. The resulting distance matrix and other results will be saved in these folders.

## Pipeline_042123.m

The `Pipeline_042123.m` script is designed to handle the complete data processing workflow. It performs the following tasks:

- Loads input data
- Pre-processes time series data
- Calculates distance matrices
- Post-processes distance matrices
- Aggregates distance matrices
- Classifies instances based on the final aggregated distance matrix

This pipeline script relies on four main functions:

- `PreProcessDat.m`: For preprocessing the data
- `CalculateDistanceMatrix.m`: For computing distance matrices
- `PostProcessDist.m`: For post-processing the distance matrices
- `AggregateDist.m`: For aggregating the distance matrices

Additionally, the pipeline script allows you to set hyperparameters, such as the window size for DTW and the number of expected clusters.




