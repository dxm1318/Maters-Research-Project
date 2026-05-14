# Maters-Research-Project

A masters research project submitted to University of Miami in partial fulfillment of the requirements for the M.S. degree in Music Engineering and Technology. The topic of the project was parameter optimization of the Synchronzied Swept Sine technique for the separation and extraction of distortion product otoacoustic emissions. This repository contains the MATLAB code used to run the simulations and a pdf copy of the paper.

The code is divided into 2 parts: preprocessing and simulation.

## Preprocessing

The preprocessing code is responsible for the generation of DPOAEs (distortion product otoacoustic emission). Using a for loop to iterate through a desginated parameter test range, 10 unique DPOAEs are generated for each parameter value. Since we are testing 3 parameters, there are 3 generation scripts: `OAE_data_r.m` , `OAE_data_alpha.m` , and `OAE_data_stim.m`. 

Each script goes through a similar process. 
  1. Prior to running the code, please ensure the amtoolbox is installed and working. This is needed to use the `verhulst2012.m` cochlear simulator.
     Run the following commands in the command window:
     `amt_start`, 
     `amt_mex`

      If there is no error, everything is working as it should. The amtoolbox is provided           in the preprocessing folder. 
 
  2. Declaration of initial parameters (sample rate, start/stop frequency for sine sweeps, etc.)
  3. Loop through parameter test range
  4. With a nested for loop, generate 10 unique DPOAEs for each parameter value using `verhulst2012.m`.   
     `SSS_gen.m` and `synchronized_swept_sine.m`  --> generate the synchronized swept sines
     `get_OAE.m` --> input synchronized swept sines into 'verhulst2012.m' to output our simulated DPOAE response     
  5. Save the preprocessed dataset as a `.mat` file: 
     `PreGenerated_OAE_R.mat`, 
     `PreGenerated_OAE_alpha.mat`, and
     `PreGenerated_OAE_Stim.mat`

## Simulation Code

The simulation code is responsible for the separation and extraction of Short-Latency and Long-Latency components for each simulated DPOAE, computing metrics, organizing and visualizing data, and statistical analysis. Each parameter has its own simulation script:
  `r_sim.m`, 
  `alpha_sim.m`, and 
  `stim_sim.m`
Similar to the preprocessing step, each of the simulation scripts follow a similar process: 
  1. If not done already, load the preprocessed data into the MATLAB workspace. In the command window type: `load('name_of_dataset.mat')`
  2. Declare initial parameters
  3. loop through parameter test range
  4. loop through the preprocessed dataset. For each element in the dataset corresponding to a single DPOAE:
      - Deconvolve the DPOAE (`VIR_deconv.m` and `synchronized_swept_sine_spectra.m`)
      - Separate the different distortion products (`synchronzied_swept_sine_IR_separation.m`)
      - for each distortion product, separate the short latency and long latnecy components(`separate_SL_LL.m`)
      - Calculate the metrics (`compute_metrics.m`): temporal separation (latency), overlap %, and Signal-to-Noise ratio of the short latency component 
  5. Organize the data into csv files: `organize_results_r.m`, `organize_results_a.m`, and `organize_results_L.m`
  6. Visualize the data (plots)
  7. Statistical Analysis (2-way ANOVA)

## Requirements to run code

  1. MATLAB 2017a or later
  2. Signal Processing Toolbox
  3. Auditory Modeling Toolbox
  4. Statistics and Machine Learning Toolbox
  5. Antonin Novak Synchronized Swept Sine Functions --> https://ant-novak.com/posts/research/2015-10-30_JAES_Swept/
  






