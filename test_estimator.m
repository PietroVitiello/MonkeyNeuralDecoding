clc
clear all
load monkeydata_training
estimator = PositionEstimator();

%%
lista = 1:300;
len = size(lista, 2);
ciak = lista(270-3:end-3)
lista(len);

%% Testing dataset creation function
[state0, eeg_train, eeg_test, x_train, x_test] = estimator.getDataset(trial, 1, 80);

%% Testing creation of dynamics matrics
[A, W, H, Q] = estimator.computeDynamics(x_train, eeg_train);

%% Testing more complex way to generate dataset
[state0, eeg_train, eeg_test, x_train, x_test] = estimator.sayonara(trial, 4, 5, 3, 80);