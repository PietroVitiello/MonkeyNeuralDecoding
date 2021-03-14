clc
clear all
load monkeydata_training
estimator = PositionEstimator();

%%
lista = 1:10;
di = diff(lista, 0, 2)

%% Testing dataset creation function
[state0, eeg_train, eeg_test, x_train, x_test] = estimator.getDataset(trial, 1, 80);

%% Testing creation of dynamics matrics
[A, W, H, Q] = estimator.computeDynamics(x_train, eeg_train);

%% Testing more complex way to generate dataset
[state0, eeg_train, eeg_test, x_train, x_test] = estimator.sayonara(trial, 4, 5, 3, 80);

%% Testing average bins
[state0, eeg_train, eeg_test, x_train, x_test] = estimator.ferromagnetico(trial, 4, 5, 3, 80);