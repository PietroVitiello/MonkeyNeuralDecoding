clc
clear all
load monkeydata_training
estimator = PositionEstimator();

%% Teting datasetcreation function
[eeg_train, eeg_test, x_train, x_test] = estimator.getLabels(trial, 1, 80);