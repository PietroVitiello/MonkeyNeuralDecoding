clc
clear all
load monkeydata_training
estimator = PositionEstimator();

%% Teting dataset creation function
[eeg_train, eeg_test, x_train, x_test] = estimator.getDataset(trial, 1, 100);

%% Testing creation of dynamics matrices
%[A, W, H, Q] = estimator.computeDynamics(x_traiestimator = PositionEstimator();n, eeg_train);

