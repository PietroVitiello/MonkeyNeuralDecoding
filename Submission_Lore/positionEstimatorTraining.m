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
  
  initial_pos = processor.get_data_matrix_initial(training_data);
  [trials, pos] = processor.get_data_matrix(training_data);
  final_pos = processor.get_data_matrix_final(training_data);
  
  start = trj.initial_positions(initial_pos);
  [avgT, stdT] = trj.averageTrajectory(pos);
  obj = trj.objective_positions(final_pos);
  
  modelParameters.initial = start;
  modelParameters.traces = avgT;
  modelParameters.deviation = stdT;
  modelParameters.objectives = obj;
  
end