function [modelParameters] = positionEstimatorTraining_8(training_data)
  % Arguments:
 
  % - training_data:
  %     training_data(n,k)              (n = trial id,  k = reaching angle)
  %     training_data(n,k).trialId      unique number of the trial
  %     training_data(n,k).spikes(i,t)  (i = neuron id, t = time)
  
  
  processor = Processing();
  a_classifier = AngleClassifier();
  
  [trials, ~] = processor.get_data_matrix(training_data);
  templates1 = a_classifier.firingTemplate(trials, 300, 1);
  par = a_classifier.neuronDistribution_mle(trials, 1, 300);
  
  modelParameters.par = par;
  modelParameters.templates1 = templates1;
  modelParameters.classifier = a_classifier;
  
end
