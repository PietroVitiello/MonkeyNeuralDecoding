function [x, y, modelParameters] = positionEstimator(test_data, modelParameters, direc)
  % Return Value:
  % - [x, y]:
  %     current position of the hand
 
  spikes = test_data.spikes;
  length_all = size(test_data.spikes,2);
  length_ = 320;
  if size(spikes, 2) <= length_
      init_spikes = sum(spikes(modelParameters.neurons, 1:length_), 2);
      angle_n = predict(modelParameters.classifier1, init_spikes');
      state0 = modelParameters.initial_params(:, angle_n);
      init_P = modelParameters.init_error_cov;
      x = test_data.startHandPos(1) + state0(1);
      y = test_data.startHandPos(2) + state0(2);
      init_x = [x; y; state0];
      modelParameters.angle_n = angle_n;
  else
      if size(spikes, 2) == 360
          init_spikes = sum(spikes(modelParameters.eccoli, 1:360), 2);
          modelParameters.angle_n = predict(modelParameters.classifier2, init_spikes');
      end
      if size(spikes, 2) == 400
          init_spikes = sum(spikes(modelParameters.eccoli, 1:400), 2);
          modelParameters.angle_n = predict(modelParameters.classifier3, init_spikes');
      end
      init_x = modelParameters.init_x;
      init_P = modelParameters.init_P;
      x = test_data.decodedHandPos(1, end);
      y = test_data.decodedHandPos(2, end);
  end
  
  A = modelParameters.A(:, :, modelParameters.angle_n);
  H = modelParameters.H(:, :, modelParameters.angle_n);
  Q = modelParameters.Q(:, :, modelParameters.angle_n);
%   inv(Q)
  det(Q);
  W = modelParameters.W(:, :, modelParameters.angle_n);
  estimator = modelParameters.pos_estimator;
  
  lag = modelParameters.lag;
  bin_size = modelParameters.bin_size;
  whos spikes
  usable_data = estimator.apply_ferromagnetico(spikes, lag, bin_size);
  
  for i = 1 : size(usable_data, 2)
      obs = usable_data(:, i);
      [init_x, init_P] = estimator.update(A, init_x, H, Q, W, init_P, obs);
      x = init_x(1); %x + (init_x(3) * bin_size);
      y = init_x(2); %y + (init_x(4) * bin_size);
  end
  
  modelParameters.init_x = init_x;
  modelParameters.init_P = init_P;
 
end