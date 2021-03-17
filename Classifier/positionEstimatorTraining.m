function [modelParameters] = positionEstimatorTraining(training_data)
  % Arguments:
 
  % - training_data:
  %     training_data(n,k)              (n = trial id,  k = reaching angle)
  %     training_data(n,k).trialId      unique number of the trial
  %     training_data(n,k).spikes(i,t)  (i = neuron id, t = time)
  
  
  processor = Processing();
  a_classifier = AngleClassifier();
 
  %silent_neurons = [8 10 11 38 49 52 73 74 76];
  %silent_neurons = [2 5 6 8 10 11 19 20 21 24 25 37 38 39 42 44 46 49 52 59 62 63 70 73 74 76 79 83 84 90 95 97];
  silent_neurons = [1 2 5 6 8 10 11 13 15 19 20 21 24 25 37 38 39 42 44 46 49 51 52 57 59 60 62 63 64 70 73 74 76 79 82 83 84 90 93 95 97];
  [clean_trial, clean_neurons] = processor.clean_dataset(training_data, silent_neurons);

  neurons_per_angle = 7;
  active_neurons = processor.mostActive(clean_trial, neurons_per_angle, clean_neurons);
  %active_neurons = processor.mostActive(training_data, neurons_per_angle, 1:98);
  eccoli_qui = 1:98; %processor.mostActive(clean_trial, 4, 300, 400);

  [samples, labels] = processor.create_dataset(training_data, active_neurons, 320, 1);
  [samples2, labels2] = processor.create_dataset(training_data, eccoli_qui, 360, 1);
  [samples3, labels3] = processor.create_dataset(training_data, eccoli_qui, 400, 1);
  
  n_neighbours = 28;
  [Mdl1, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples, labels);
  [Mdl2, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples2, labels2);
  [Mdl3, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples3, labels3);
  
  modelParameters.classifier1 = Mdl1;
  modelParameters.classifier2 = Mdl2;
  modelParameters.classifier3 = Mdl3;
  modelParameters.neurons = active_neurons;
  modelParameters.eccoli = eccoli_qui;
  
end
