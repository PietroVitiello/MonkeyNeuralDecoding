clc
clear all
load monkeydata_training
processor = Processing();
a_classifier = AngleClassifier();

%% KNN
silent_neuron = [8 10 11 38 49 52 73 74 76];
clean_trial = processor.clean_dataset(trial, silent_neuron);

runs = 100;
neighbours = 1:5:40;
n_correct_angles = zeros(length(neighbours), runs);
angle_distribution=zeros(runs, size(clean_trial,2));

figure()
hold on;
for i=1:runs
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
    
    for k = 1:length(neighbours)
        [train_metrics, test_metrics] = a_classifier.knn_classifier(training_samples, training_labels, test_samples, test_labels, k);
        n_correct_angles(k,i) = test_metrics;
        plot(n_correct_angles(k, :))
    end
%     classified_data = a_classifier.k_nn(training_samples, training_labels, test_samples, 28);
%     correct_angles = classified_data == test_labels.';
%     n_correct_angles(i) = sum(correct_angles); %Number of correctly classified trials

end

correct_angles = classified_data == test_labels.';
n_correct_angles = sum(correct_angles); %Number of correctly classified trials

%% MLE
silent_neuron = [8 10 11 38 49 52 73 74 76];
clean_trial = processor.clean_dataset(trial, silent_neuron);
active_neurons = processor.mostActive(clean_trial, 4);

[train_mx, test_mx] = processor.data_as_matrix(clean_trial, active_neurons, 80);

[estimated_angles, true_angles] = a_classifier.likelihood(train_mx, test_mx);

correct_angles = estimated_angles == true_angles;
accuracy = sum(correct_angles)/length(true_angles)*100

%% Covariance
silent_neuron = [8 10 11 38 49 52 73 74 76];
clean_trial = processor.clean_dataset(trial, silent_neuron);
active_neurons = processor.mostActive(clean_trial, 7);

[train_mx, ~] = processor.data_as_matrix(clean_trial, active_neurons, 80);

covariance_matrix = processor.covariance(train_mx, 1)

%% Multidimensional MLE
silent_neuron = [8 10 11 38 49 52 73 74 76];
clean_trial = processor.clean_dataset(trial, silent_neuron);
active_neurons = processor.mostActive(clean_trial, 9);

[train_mx, test_mx] = processor.data_as_matrix(clean_trial, active_neurons, 90);

[estimated_angles, true_angles] = a_classifier.multidimensional_mle(train_mx, test_mx);

correct_angles = estimated_angles == true_angles;
accuracy = sum(correct_angles)/length(true_angles)*100
