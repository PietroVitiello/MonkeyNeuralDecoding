clc
clear all
load monkeydata_training
%% 
silent_neuron = [8 10 11 38 49 52 73 74 76];
processor = Processing();
clean_trial = processor.clean_dataset(trial, silent_neuron);

runs = 100;
n_correct_angles = zeros(1, runs);
angle_distribution=zeros(runs, size(clean_trial,2));
for i=1:runs
    i
    active_neurons = processor.mostActive(clean_trial, 9);
    [samples, labels] = processor.create_dataset(trial, active_neurons);
    
    training_samples = samples(81:800, :);
    training_labels = labels(81:800, :);
    test_samples = samples(1:80, :);
    test_labels = labels(1:80, :);
    
    for angle=1:8
        idx = [];
        idxs = find(test_labels==angle);
        angle_distribution(i, angle) = (length(idxs))/(size(test_labels,1));
    end

    a_classifier = AngleClassifier();
    neighbours = 1:40;

    [train_metrics, test_metrics] = a_classifier.knn_classifier(training_samples, training_labels, test_samples, test_labels, 28);
    n_correct_angles(i) = test_metrics;
    
%     classified_data = a_classifier.k_nn(training_samples, training_labels, test_samples, 28);
%     correct_angles = classified_data == test_labels.';
%     n_correct_angles(i) = sum(correct_angles); %Number of correctly classified trials

end

plot(n_correct_angles)




