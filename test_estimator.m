clc
clear all
load monkeydata_training

estimator = PositionEstimator();

%% Testing dataset creation function
[state0, eeg_train, eeg_test, x_train, x_test] = estimator.getDataset(trial, 1, 80);
size(x_train{1})

%% Testing creation of dynamics matrics
[A, W, H, Q] = estimator.computeDynamics(x_train, eeg_train);
