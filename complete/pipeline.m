%% Nature Communications Paper Analysis Pipeline
% Darren Price 2016

% Epoch the data 
edit copy_epoch_data_avpassive.m
edit copy_epoch_data_avtask.m

% Calculate principal components
edit PCA_Sensors_allsubs_AVPASSIVE.m % Figure 1 and Supp Fig 2
edit PCA_Sensors_allsubs_AVPASSIVE_0_150ms.m % Supp Fig 1 (reruns fitting with short time windows)
edit PCA_Sensors_allsubs_AVTASK.m % Supp Fig 2
edit PCA_Sensors_allsubs_AVTASK_ApplyPassiveWeights.m

% Plot Delay Estimation Figures
edit Stats_ERPs_Paper.m
edit Stats_ERPs_Paper_shortwin.m % Supplementary Figure 1

% Run mediation analysis
edit run_robust_mediation.m

% Answers to Reviewer's Comments 
edit runtest_gaussian_peaks_2Dgrid.m % Supp Figure 3
edit runtest_real_data_peak_vs_template_avpassive.m % Supp Figure 4
edit create_powercalc_simulations.m % Supp Figure 5
