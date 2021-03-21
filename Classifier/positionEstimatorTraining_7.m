function [modelParameters] = positionEstimatorTraining_7(training_data)

  processor = Processing();
  a_classifier = AngleClassifier();
  
  neurons_per_angle = 12;
  [~, ~] = processor.mostActive(training_data, neurons_per_angle, 1:98);
  active_neurons = 1:98;
  [samples, labels] = processor.create_dataset(training_data, active_neurons, 320, 1); % NOT USING "active_neurons"
  
  [trials, ~] = processor.get_data_matrix(training_data);
  templates1 = a_classifier.firingTemplate(trials, 320, 1);
  av_diffs = processor.get_dataset_2(templates1, trials, active_neurons, 320, 1);
  samples = [samples av_diffs.*20];
  
  temp = [samples labels];
  rand_d = temp(randperm(size(temp, 1)), :);
  samples = rand_d(:, 1:size(samples, 2));
  labels = rand_d(:, size(samples, 2)+1);
  
  size(samples)
  
  n_neighbours = 28;
  Mdl = a_classifier.knn_classifier(n_neighbours, samples, labels);
  
  modelParameters.templates1 = templates1;
  modelParameters.classifier = a_classifier;
  modelParameters.neurons = active_neurons;
  modelParameters.knn = Mdl;

end
