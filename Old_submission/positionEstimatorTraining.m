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
  [samples, labels] = processor.create_dataset(training_data, active_neurons);
  
  n_neighbours = 28;
  [Mdl, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples, labels);
  
  lag = 1;
  bin_size = 3;
  [state0, eeg_train, ~, x_train, ~] = estimator.sayonara(training_data, lag, bin_size, 2, 100);
  
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
  modelParameters.classifier = Mdl;
  modelParameters.initial_params = state0;
  modelParameters.lag = lag;
  modelParameters.bin_size = bin_size;
  modelParameters.neurons = active_neurons;
  modelParameters.pos_estimator = estimator;
  modelParameters.init_error_cov = ones(6,6);
end
