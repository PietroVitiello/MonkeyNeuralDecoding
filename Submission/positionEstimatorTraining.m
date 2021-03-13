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
  %     - neurons
  %     - pos_estimator
  %     - init_error_cov
  
  processor = Processing();
  a_classifier = AngleClassifier();
  estimator = PositionEstimator_cl();
 
  silent_neurons = [8 10 11 38 49 52 73 74 76]; 
  %9 neurons with the least firing across trials and angles
  clean_trial = processor.clean_dataset(training_data, silent_neurons);
  %Silent neurons are cleared out of the training data
  
  neurons_per_angle = 9;
  neurons_classifier = processor.mostActive(clean_trial, neurons_per_angle); 
  %Most active neurons across trials and angles
  [samples, labels] = processor.create_dataset(training_data, neurons_classifier);
  %Sample & label dataset is created to build a K-nearest neighbours
  %classifier
  
  n_neighbours = 28; %Number of nearest neighbours
  [Mdl, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples, labels);
  %Mdl contains the parameters of the classifier
  
  neurons_per_angle = 3;
  neurons_estimator = processor.mostActive(clean_trial, neurons_per_angle, 300, 571); 
 
  [state0, eeg_train, ~, x_train, ~] = estimator.getDataset(training_data, 1, 100);
  %state0 is the average velocity of the hand during the 300ms prior to the
  %movement; eeg_train is the training dataset of spike trains from 300 ms
  %onwards; x_train is the training dataset of hand velocity from 300 ms
  %onwards
  
  eeg_train = estimator.non_redundant(eeg_train, neurons_estimator);
   
  [A, W, H, Q] = estimator.computeDynamics(x_train, eeg_train);
  %A, W, H, Q are the Kalman filter matrices
  
  modelParameters.A = A; 
  modelParameters.W = W;
  modelParameters.H = H;
  modelParameters.Q = Q;
  modelParameters.classifier = Mdl;
  modelParameters.initial_params = state0;
  modelParameters.neurons_classifier = neurons_classifier;
  modelParameters.neurons_estimator = neurons_estimator;
  modelParameters.pos_estimator = estimator;
  modelParameters.init_error_cov = ones(2,2);
end
