function [modelParameters] = positionEstimatorTraining(training_data)
processor = Processing();
a_classifier = AngleClassifier();
estimator = PositionEstimator_cl();

%     silent_neurons = [8 10 11 38 49 52 73 74 76];
%     clean_trial = processor.clean_dataset(training_data, silent_neurons);

neurons_per_angle = 9;
active_neurons = processor.mostActive(training_data, neurons_per_angle);

[samples, labels] = processor.create_dataset(training_data, active_neurons, 320, 1);
[samples2, labels2] = processor.create_dataset(training_data, 1:98, 360, 1);
[samples3, labels3] = processor.create_dataset(training_data, 1:98, 400, 1);

n_neighbours = 28;
[Mdl1, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples, labels);
[Mdl2, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples2, labels2);
[Mdl3, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples3, labels3);

lag = 5;
bin_size = 5;
order = 3;
[state0, eeg_train, ~, x_train, ~] = estimator.ferromagnetico(training_data, lag, bin_size, order, 100);

[A, W, H, Q] = estimator.computeDynamics(x_train, eeg_train);

modelParameters.A = A; 
modelParameters.W = W;
modelParameters.H = H;
modelParameters.Q = Q;
modelParameters.classifier1 = Mdl1;
modelParameters.classifier2 = Mdl2;
modelParameters.classifier3 = Mdl3;
modelParameters.initial_params = state0;
modelParameters.lag = lag;
modelParameters.bin_size = bin_size;
modelParameters.neurons = active_neurons;
modelParameters.pos_estimator = estimator;
modelParameters.init_error_cov = ones((order+1)*2);

end
