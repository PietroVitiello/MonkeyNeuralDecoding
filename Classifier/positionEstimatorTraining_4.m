function [modelParameters] = positionEstimatorTraining_4(training_data)

  processor = Processing();
  a_classifier = AngleClassifier();
  [trials, ~] = processor.get_data_matrix(training_data);
  
  neuronsXangle = 8;
  [pref_mag, pref_neuron, sum_activity] = processor.anglePreference(trials, 1, 300);
  templates = a_classifier.firingTemplate(trials, 300, 1);
  angle_distributions = processor.firingDistribution(trials, 1, 320);
  pref_mag = pref_mag(1:neuronsXangle, :);
  pref_neuron = pref_neuron(1:neuronsXangle, :);
  
  modelParameters.templates = templates;
  modelParameters.distributions = angle_distributions;
  modelParameters.pref_mag = pref_mag;
  modelParameters.pref_neuron = pref_neuron;
  modelParameters.sum_activity = sum_activity;
  modelParameters.classifier = a_classifier;
  
end
