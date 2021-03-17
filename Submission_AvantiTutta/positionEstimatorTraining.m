function [modelParameters] = positionEstimatorTraining(training_data)
  % Arguments:
 
  % - training_data:
  %     training_data(n,k)              (n = trial id,  k = reaching angle)
  %     training_data(n,k).trialId      unique number of the trial
  %     training_data(n,k).spikes(i,t)  (i = neuron id, t = time)
  %     training_data(n,k).handPos(d,t) (d = dimension [1-3], t = time)
  
  % Return Value:
  
  % - modelParameters:
  %     - A
  %     - W
  %     - H
  %     - Q
  %     - classifier
  %     - initial parameters
  
  
  processor = Processing();
  a_classifier = AngleClassifier();
  estimator = PositionEstimator_cl();
 
  silent_neurons = [8 10 11 38 49 52 73 74 76];
  clean_trial = processor.clean_dataset(training_data, silent_neurons);
  
  neurons_per_angle = 9;
  active_neurons = processor.mostActive(clean_trial, neurons_per_angle);
  eccoli_qui = 1:98; processor.mostActive(clean_trial, 4, 300, 400);
  
  [samples, labels] = processor.create_dataset(training_data, active_neurons, 320, 1);
  [samples2, labels2] = processor.create_dataset(training_data, eccoli_qui, 360, 1);
  [samples3, labels3] = processor.create_dataset(training_data, eccoli_qui, 400, 1);
  
  n_neighbours = 28;
  [Mdl1, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples, labels);
  [Mdl2, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples2, labels2);
  [Mdl3, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples3, labels3);
  
  lag = 5;
  bin_size = 5;
  order = 3;
  [state0, eeg_train, ~, x_train, ~] = estimator.ferromagnetico_2(training_data, lag, bin_size, order, 100);
  
  assignin('base', 'eeg_train', eeg_train);
   
  [A, W, H, Q] = estimator.computeDynamics(x_train, eeg_train);
  assignin('base', 'A', A);
  assignin('base', 'W', W);
  assignin('base', 'H', H);
  assignin('base', 'Q', Q);
  
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
  modelParameters.eccoli = eccoli_qui;
  modelParameters.pos_estimator = estimator;
  modelParameters.init_error_cov = ones((order+1)*2);
end
