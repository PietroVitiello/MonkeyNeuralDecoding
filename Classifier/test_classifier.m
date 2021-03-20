clc
clear all
load monkeydata0
classifier = AngleClassifier();
process = Processing();

% Create matrices
[trials, pos] = process.get_data_matrix(trial);

%% Template of average activity
templates = classifier.firingTemplate(trials, 300, 1)';

%% Find most specific neurons
specific_neurons = process.mostSpecific(trials, 1, 320, 98);

%% Find most specific neurons
distributions = process.firingDistribution(trials, 1, 320);

%% Find preference
[pref_mag, pref_neuron, ~] = process.anglePreference(trials, 1, 320);

%% test neuron distribution mle
par = classifier.neuronDistribution_mle(trials, 1, 300);

%% uuu
A = reshape(1:125, 5,5,5);
% i1 = repmat(1:3,3,1)
% i2 = repmat([1:3]',1,3)
i1 = 1:2
i2 = 1:3
B = squeeze(A(i2, i1, 1))

