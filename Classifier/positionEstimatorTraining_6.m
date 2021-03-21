function [modelParameters] = positionEstimatorTraining_6(training_data)
  % Arguments:
 
  % - training_data:
  %     training_data(n,k)              (n = trial id,  k = reaching angle)
  %     training_data(n,k).trialId      unique number of the trial
  %     training_data(n,k).spikes(i,t)  (i = neuron id, t = time)
  
  
  processor = Processing();
  a_classifier = AngleClassifier();
 
  %silent_neurons = [8 10 24 25 28 38 42 46 49 52 73 74 76 79 83]; %Threshold of 5
  %silent_neurons = [8 10 11 19 20 24 25 28 38 39 42 46 49 52 54 58 62 73 74 76 79 83 84 95]; %Threshold of 10
  %silent_neurons = [2 5 6 8 10 11 19 20 21 23 24 25 26 28 30 35 37 38 39 42 44 46 49 52 53 54 58 59 60 62 63 64 70 73 74 76 79 82 83 84 90 95 97]; %Threshold of 15
  %[clean_trial, clean_neurons] = processor.clean_dataset(training_data, silent_neurons);

  neurons_per_angle = 10;
  %[active_neurons active_neurons_matrix] = processor.mostActive(clean_trial, neurons_per_angle, clean_neurons);
  [active_neurons, active_neurons_matrix] = processor.mostActive(training_data, neurons_per_angle, 1:98);
  eccoli_qui = 1:98; %processor.mostActive(clean_trial, 4, 300, 400);

  [samples, labels] = processor.create_dataset(training_data, eccoli_qui, 320, 1);
  [samples2, labels2] = processor.create_dataset(training_data, eccoli_qui, 360, 1);
  [samples3, labels3] = processor.create_dataset(training_data, eccoli_qui, 400, 1);
  
  n_neighbours = 28;
  [Mdl1, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples, labels);
  [Mdl2, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples2, labels2);
  [Mdl3, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples3, labels3);
  
  [trials, ~] = processor.get_data_matrix(training_data);
  templates = a_classifier.firingTemplate(trials, 300, 1);
  
  av_diffs = processor.get_dataset_2(templates, trials, eccoli_qui, 320, 1);
  samples = [samples av_diffs.*20];
  
  temp = [samples labels];
  rand_d = temp(randperm(size(temp, 1)), :);
  samples = rand_d(:, 1:size(samples, 2));
  labels = rand_d(:, size(samples, 2)+1);
  [Mdl4, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples, labels);
  
  angle_distributions = processor.firingDistribution(trials, 1, 320);
  
  par = a_classifier.neuronDistribution_mle(trials, 1, 300);
  
  
  average_spike_trains = processor.averageTrial(trials, active_neurons, 1, 320);
  
  modelParameters.classifier1 = Mdl1;
  modelParameters.classifier2 = Mdl2;
  modelParameters.classifier3 = Mdl3;
  modelParameters.classifier4 = Mdl4;
  modelParameters.neurons = active_neurons;
  modelParameters.neuron_matrix = active_neurons_matrix;
  modelParameters.classifier = a_classifier;
  modelParameters.eccoli = eccoli_qui;
  modelParameters.templates = templates;
  modelParameters.distributions = angle_distributions;
  modelParameters.par = par;
  modelParameters.training_average = average_spike_trains;
  
end
