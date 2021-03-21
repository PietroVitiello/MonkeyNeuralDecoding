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

%% Test passing vectors as index
A = reshape(1:125, 5,5,5);
% i1 = repmat(1:3,3,1)
% i2 = repmat([1:3]',1,3)
i1 = 1:2;
i2 = 1:3;
B = squeeze(A(i2, i1, 1));

%% multivariate gaussian mle
silent_neuron = [8 24 25 38 42 49 52 54 73 74 76];
% silent_neuron = [8 10 11 38 49 52 73 74 76];
clean_trial = process.clean_dataset(trial, []);
% active_neurons = process.mostActive(clean_trial, 9);

[train_mx, test_mx] = process.data_as_matrix(clean_trial, 1:98, 90);
train_mx(:,:,silent_neuron,:) = [];
test_mx(:,:,silent_neuron,:) = [];
train_mx(:,:,1:30,:) = [];
test_mx(:,:,1:30,:) = [];

[estimated_angles, true_angles] = classifier.multidimensional_mle(train_mx, test_mx);

correct_angles = estimated_angles == true_angles;
accuracy = sum(correct_angles)/length(true_angles)*100

%% Testing reshaping
A = repmat([1:3]', 1, 3);
A = cat(3, A, repmat([4:6]', 1, 3));
B = reshape(permute(A, [2 1 3]), 9, 2);

A = permute(reshape([1:16], 2, 4, 2), [2 1 3])
B = permute(reshape(permute(A, [2 1 3]), 2, 4, 3), [2 1 3]);

%% Testing mulivariate gaussian version 2
% silent_neuron = [8 24 25 38 42 49 52 54 73 74 76];
silent_neuron = [8 24];
clean_trial = process.clean_dataset(trial, []);
% active_neurons = process.mostActive(clean_trial, 9);

[train_mx, test_mx] = process.data_as_matrix(clean_trial, 1:98, 90);
% train_mx(:,:,silent_neuron,:) = [];
% test_mx(:,:,silent_neuron,:) = [];

[estimated_angles, true_angles] = classifier.multidimensional_mlev2(train_mx, test_mx, 2);

correct_angles = estimated_angles == true_angles;
accuracy = sum(correct_angles)/length(true_angles)*100





