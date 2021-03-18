clc
clear all
load monkeydata0
classifier = AngleClassifier();
process = Processing();

% Create matrices
[trials, pos] = process.get_data_matrix(trial);

%% Find most specific neurons
specific_neurons = process.mostSpecific(trials, 1, 320, 98);

%% Find most specific neurons
distributions = process.firingDistribution(trials, 1, 320);


