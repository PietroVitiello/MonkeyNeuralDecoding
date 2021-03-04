clc
clear all
load monkeydata_training
%% 
silent_neuron = [8 10 11 38 49 52 73 74 76];
processor = Processing();
clean_trial = processor.clean_dataset(trial, silent_neuron);

active_neurons = processor.mostActive(clean_trial, 4);

[samples, labels] = processor.create_dataset(trial, active_neurons);

training_samples = samples(81:800, :);
training_labels = labels(81:800, :);
test_samples = samples(1:80, :);
test_labels = labels(1:80, :);

a_classifier = AngleClassifier();
classified_data = a_classifier.k_nn(training_samples, training_labels, test_samples, 28);

correct_angles = classified_data == test_labels.';
n_correct_angles = sum(correct_angles); %Number of correctly classified trials