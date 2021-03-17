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
  estimator = PositionEstimator_cl();
  
  lag = 5;
  bin_size = 5;
  order = 3;
  [state0, eeg_train, ~, x_train, ~] = estimator.ferromagnetico_2(training_data, lag, bin_size, order, 100);
  
  assignin('base', 'eeg_train', eeg_train);
   
  [A, W, H, Q] = estimator.computeDynamics(x_train, eeg_train);
  assignin('base', 'A', A);
  assignin('base', 'W', W);
  assignin('base', 'H', H);
  assignin('base', 'Q', Q);
  
  modelParameters.A = A; 
  modelParameters.W = W;
  modelParameters.H = H;
  modelParameters.Q = Q;
  modelParameters.initial_params = state0;
  modelParameters.lag = lag;
  modelParameters.bin_size = bin_size;
  modelParameters.pos_estimator = estimator;
  modelParameters.init_error_cov = ones((order+1)*2);
end
