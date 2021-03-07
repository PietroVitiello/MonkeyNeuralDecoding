clc
clear all
load monkeydata_training

estimator = PositionEstimator();

%% Teting dataset creation function
[eeg_train, eeg_test, x_train, x_test] = estimator.getLabels(trial, 1, 80);

%% Testing creation of dynamics matrices
[A, W, H, Q] = estimator.computeDynamics(x_train, eeg_train);
