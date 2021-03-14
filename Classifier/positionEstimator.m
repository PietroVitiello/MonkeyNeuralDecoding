function [angle, modelParameters] = positionEstimator(test_data, modelParameters)
 
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
  end
 
end