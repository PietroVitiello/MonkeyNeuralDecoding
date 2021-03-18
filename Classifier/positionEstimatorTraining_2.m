function [modelParameters] = positionEstimatorTraining_2(training_data)
  % Arguments:
 
  % - training_data:
  %     training_data(n,k)              (n = trial id,  k = reaching angle)
  %     training_data(n,k).trialId      unique number of the trial
  %     training_data(n,k).spikes(i,t)  (i = neuron id, t = time)
  
  
  processor = Processing();
  a_classifier = AngleClassifier();
  
  [trials, ~] = processor.get_data_matrix(training_data);
  templates1 = a_classifier.firingTemplate(trials, 300, 1);
%   templates2 = a_classifier.firingTemplate(trials, 400, 1);
  
  modelParameters.templates1 = templates1;
%   modelParameters.templates2 = templates2;
  modelParameters.classifier = a_classifier;
  
end
