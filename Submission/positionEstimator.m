function [x, y, modelParameters] = positionEstimator(test_data, modelParameters)
  % Return Value:
  % - [x, y]:
  %     current position of the hand
 
  spikes = test_data.spikes;
  length_all = size(test_data.spikes,2);
  length_ = 320;
  if size(spikes, 2) <= length_
      init_spikes = sum(spikes(modelParameters.neurons, 1:length_), 2);
      angle_n = predict(modelParameters.classifier, init_spikes');
      init_x = modelParameters.initial_params;
      init_P = modelParameters.init_error_cov;
      x = test_data.startHandPos(1);
      y = test_data.startHandPos(2);
      modelParameters.angle_n = angle_n;
  else
      init_x = modelParameters.init_x;
      init_P = modelParameters.init_P;
      x = test_data.decodedHandPos(1, end);
      y = test_data.decodedHandPos(2, end);
  end
  
  A = modelParameters.A(:, :, modelParameters.angle_n);
  H = modelParameters.H(:, :, modelParameters.angle_n);
  Q = modelParameters.Q(:, :, modelParameters.angle_n);
  inv(Q);
  W = modelParameters.W(:, :, modelParameters.angle_n);
  estimator = modelParameters.pos_estimator;
  
  for i = length_all-19
      obs = spikes(:, i);
      [init_x, init_P] = estimator.update(A, init_x, H, Q, W, init_P, obs);
      x = x + init_x(1);
      y = y + init_x(2);
  end
  
  modelParameters.init_x = init_x;
  modelParameters.init_P = init_P;
 
end