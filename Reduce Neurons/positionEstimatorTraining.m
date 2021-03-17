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
  
  neurons_per_angle = 8;
  neurons_classifier = processor.mostActive(clean_trial, neurons_per_angle);
  %neurons_classifier = processor.mostActive(training_data, neurons_per_angle);
  [samples, labels] = processor.create_dataset(training_data, neurons_classifier, 320, 1);
  [samples2, labels2] = processor.create_dataset(training_data, [1:98], 360, 1);
  [samples3, labels3] = processor.create_dataset(training_data, [1:98], 400, 1);
  
  n_neighbours = 28;
  [Mdl1, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples, labels);
  [Mdl2, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples2, labels2);
  [Mdl3, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples3, labels3);
  
  neurons_per_angle = 80;
  mode = 'all';
  if strcmp(mode, 'matrix')
      neurons_estimator_matrix = processor.mostActive(training_data, neurons_per_angle, 300, 571, mode);
  end
  if strcmp(mode, 'vector')
      neurons_estimator = processor.mostActive(training_data, neurons_per_angle, 300, 571, mode); 
      neurons_estimator_matrix = zeros(length(neurons_estimator), size(training_data, 2));
      for angle_n = 1:size(training_data, 2)
          neurons_estimator_matrix(:, angle_n) = neurons_estimator';
      end
  end
  if strcmp(mode, 'all')
      neurons = [1 2 3 4 5 6 7 9 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 39 40 41 42 43 44 45 46 47 48 50 51 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 75 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98];
      
%       neurons_estimator_matrix = zeros(size(training_data(1, 1).spikes, 1), size(training_data, 2));
%       for column = 1:size(neurons_estimator_matrix)
%           neurons_estimator_matrix(:, column) = 1:98;
%       end
      neurons_estimator_matrix = zeros(size(neurons, 2), size(training_data, 2));
      for column = 1:size(neurons_estimator_matrix)
          neurons_estimator_matrix(:, column) = neurons';
      end
  end
  neurons_estimator_matrix;
  
  lag = 5;
  bin_size = 5;
  [state0, eeg_train, ~, x_train, ~] = estimator.ferromagnetico(training_data, lag, bin_size, 3, 100);
  
  eeg_train = estimator.non_redundant(eeg_train, neurons_estimator_matrix);
  eeg_train{1, 1}
  
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
  modelParameters.neurons_classifier = neurons_classifier;
  modelParameters.neurons_estimator = neurons_estimator_matrix;
  modelParameters.pos_estimator = estimator;
  modelParameters.init_error_cov = ones(8,8);
end
