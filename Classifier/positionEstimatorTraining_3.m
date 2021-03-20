function [modelParameters] = positionEstimatorTraining_3(training_data)
  % Arguments:
 
  % - training_data:
  %     training_data(n,k)              (n = trial id,  k = reaching angle)
  %     training_data(n,k).trialId      unique number of the trial
  %     training_data(n,k).spikes(i,t)  (i = neuron id, t = time)
  
  
  processor = Processing();
  a_classifier = AngleClassifier();
  
  [trials, ~] = processor.get_data_matrix(training_data);
  templates = a_classifier.firingTemplate(trials, 300, 1);
  angle_distributions = processor.firingDistribution(trials, 1, 320);
  
  modelParameters.classifier = a_classifier;
  modelParameters.templates = templates;
  modelParameters.distributions = angle_distributions;
  
end
