function [modelParameters] = positionEstimatorTraining(training_data)
  % Arguments:
  
  % - training_data:
  %     training_data(n,k)              (n = trial id,  k = reaching angle)
  %     training_data(n,k).trialId      unique number of the trial
  %     training_data(n,k).spikes(i,t)  (i = neuron id, t = time)
  %     training_data(n,k).handPos(d,t) (d = dimension [1-3], t = time)
  
  % ... train your model
  
  % Return Value:
  
  % - modelParameters:
  %     single structure containing all the learned parameters of your
  %     model and which can be used by the "positionEstimator" function.
  
  processor = Processing();
  a_classifier = AngleClassifier();
  estimator = PositionEstimator();
  
  silent_neuron = [8 10 11 38 49 52 73 74 76];
  clean_trial = processor.clean_dataset(trial, silent_neuron);
  
  neurons_per_angle = 9;
  active_neurons = processor.mostActive(clean_trial, neurons_per_angle);
  [samples, labels] = processor.create_dataset(trial, active_neurons);
  
  n_neighbours = 28;
  [Mdl, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples, labels);
  
  [state0, eeg_train, ~, x_train, ~] = estimator.getDataset(trial, 1, 80);
   
  [A, W, H, Q] = estimator.computeDynamics(x_train, eeg_train);
  
  modelParameters = [A, W, H, Q, Mdl, state0];
  
end
