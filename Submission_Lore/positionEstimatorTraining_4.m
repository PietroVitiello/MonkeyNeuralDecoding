function [modelParameters] = positionEstimatorTraining_4(training_data)
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
  
  initial_pos = processor.get_data_matrix_initial(training_data);
  [trials, pos] = processor.get_data_matrix(training_data);
  final_pos = processor.get_data_matrix_final(training_data);
  
  templates1 = a_classifier.firingTemplate_3D(trials, 300, 1, 150);
  
%   templates = a_classifier.firingTemplate(trials, 300, 1);
%   angle_distributions = processor.firingDistribution(trials, 1, 320);
%   
%   [samples, labels] = processor.create_dataset(training_data, 1:98, 320, 1);
%   [samples2, labels2] = processor.create_dataset(training_data, 1:98, 360, 1);
%   [samples3, labels3] = processor.create_dataset(training_data, 1:98, 400, 1);
%   n_neighbours = 28;
%   [Mdl1, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples, labels);  
%   [Mdl2, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples2, labels2);
%   [Mdl3, ~, ~] = a_classifier.knn_classifier(n_neighbours, samples3, labels3);
  
  start = trj.initial_positions(initial_pos);
  [avgT, stdT] = trj.averageTrajectory(pos);
  obj = trj.objective_positions(final_pos);
  
  modelParameters.initial = start;
  modelParameters.traces = avgT;
  modelParameters.deviation = stdT;
  modelParameters.objectives = obj;
  modelParameters.classifier = a_classifier;
%   modelParameters.classifier1 = Mdl1;
%   modelParameters.classifier2 = Mdl2;
%   modelParameters.classifier3 = Mdl3;
%   modelParameters.templates = templates;
%   modelParameters.distributions = angle_distributions;
  modelParameters.templates1 = templates1;
  
end