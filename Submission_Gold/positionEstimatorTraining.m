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
  trj = Trajectory();
  a_classifier = AngleClassifier();
  
  [trials, ~] = processor.get_data_matrix(training_data);
  templates = a_classifier.firingTemplate(trials, 300, 1);
   
  [~, pos] = processor.get_data_matrix(training_data);
  obj = trj.objective_positions(pos);
  [avgT, ~] = trj.averageTrajectory(pos);
  
  modelParameters.templates = templates;
  modelParameters.classifier = a_classifier;
  modelParameters.objectives = obj;
  modelParameters.traces = avgT;
  
end
