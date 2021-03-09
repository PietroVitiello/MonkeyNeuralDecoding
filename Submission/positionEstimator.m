function [x, y, modelParameters] = positionEstimator(test_data, modelParameters)
  % Return Value:
  
  % - [x, y]:
  %     current position of the hand
  
  spikes = test_data.spikes;
  [x,y] = test_data.decodedHandPos;
  length_all = size(test_data.spikes,2);
  length_ = 320;
  if size(spikes, 2) <= length_
      init_spikes = spikes(modelParameters.neurons, 1:length_);
      angle_n = predict(Mdl, init_spikes);
      init_x = modelParameters.initial_params;
      init_P = modelParameters.init_P;
      [x,y] = test_data.startHandPos;
  end
  A = modelParameters.A(:, :, angle_n);
  H = modelParameters.H(:, :, angle_n);
  Q = modelParameters.Q(:, :, angle_n);
  W = modelParameters.W(:, :, angle_n);
  estimator = modelParameters.pos_estimator;
  
  for i = length_all-19:length_all
      obs = spikes(:, i);
      [init_x, init_P] = estimator.update(A, init_x, H, Q, W, init_P, obs);
      x = x + init_x(1);
      y = y + init_x(2);
  end
 
end