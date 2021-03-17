function [x, y, modelParameters] = positionEstimator_v2(test_data, modelParameters)

      spikes = test_data.spikes;
      length_ = 320;
      redo = 0;
      if size(spikes, 2) <= length_
          init_spikes = sum(spikes(modelParameters.neurons, 1:length_), 2);
          angle_n = predict(modelParameters.classifier1, init_spikes');
          state0 = modelParameters.initial_params(:, angle_n);
          init_P = modelParameters.init_error_cov;
          modelParameters.startingHand = [test_data.startHandPos(1), test_data.startHandPos(2)];
          x = test_data.startHandPos(1) + state0(1);
          y = test_data.startHandPos(2) + state0(2);
          init_x = [x; y; state0];
          modelParameters.angle_n = angle_n;
      else
          if size(spikes, 2) == 360
              init_spikes = sum(spikes(1:98, 1:360), 2);
              temp_angle = predict(modelParameters.classifier2, init_spikes');
%               modelParameters.angle_n = predict(modelParameters.classifier2, init_spikes');
              if (temp_angle ~= modelParameters.angle_n)
                  redo = 1;
                  modelParameters.angle_n = temp_angle;
              end
          elseif size(spikes, 2) == 400
              init_spikes = sum(spikes(1:98, 1:400), 2);
              temp_angle = predict(modelParameters.classifier3, init_spikes');
%               modelParameters.angle_n = predict(modelParameters.classifier3, init_spikes');
              if (temp_angle ~= modelParameters.angle_n)
                  redo = 1;
                  modelParameters.angle_n = temp_angle;
              end
          end
          if redo == 0
              init_x = modelParameters.init_x;
              init_P = modelParameters.init_P;
              x = test_data.decodedHandPos(1, end);
              y = test_data.decodedHandPos(2, end);
          else
              x = modelParameters.startingHand(1);
              y = modelParameters.startingHand(2);
              state0 = modelParameters.initial_params(:, temp_angle);
              init_x = [x+state0(1); y+state0(2); state0];
              init_P = modelParameters.init_error_cov;
          end
      end

      A = modelParameters.A(:, :, modelParameters.angle_n);
      H = modelParameters.H(:, :, modelParameters.angle_n);
      Q = modelParameters.Q(:, :, modelParameters.angle_n);
      W = modelParameters.W(:, :, modelParameters.angle_n);
      estimator = modelParameters.pos_estimator;

      lag = modelParameters.lag;
      bin_size = modelParameters.bin_size;
%       usable_data = estimator.apply_ferromagnetico_temp(spikes, lag, bin_size, redo);
      usable_data = estimator.apply_ferromagnetico(spikes, lag, bin_size);
      
      for i = 1 : size(usable_data, 2)
          obs = usable_data(:, i);
          [init_x, init_P] = estimator.update(A, init_x, H, Q, W, init_P, obs);
          x = init_x(1);
          y = init_x(2);
      end

      modelParameters.init_x = init_x;
      modelParameters.init_P = init_P;
 
end