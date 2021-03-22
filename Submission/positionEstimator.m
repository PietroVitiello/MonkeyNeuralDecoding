function [x, y, modelParameters] = positionEstimator(test_data, modelParameters)

  
  len = size(test_data.spikes,2);
  time_point = (len-300)/20;
  noise = zeros(2, 1);
  scale = zeros(2, 1);
  
  if time_point == 1
     
     templates = modelParameters.templates2;
     spikes = test_data.spikes;
     
     n_n = size(spikes, 1);
     n_bin = 300/100;
     spikes = spikes(:,1:300);

     B = movsum(spikes, 100, 2, 'Endpoints','discard');
     B = reshape(B(:,1:100:end), n_n*n_bin, 1);
     B = [mean(spikes, 2) ; B];
     B = [B(1:n_n) ./ repmat(sum( ...
        B(1:n_n)), n_n, 1); B];

     [~, angle] = min(sum(abs(templates - ...
                 repmat(B, 1, 8)' ...
                 ), 2));
      
      modelParameters.angle_n = mode(angle);
      
      if test_data.startHandPos(1, 1) < modelParameters.initial(1, modelParameters.angle_n)
          noise(1, 1) = -1;
          scale(1, 1) = abs(test_data.startHandPos(1, 1) - modelParameters.initial(1, modelParameters.angle_n));
      elseif test_data.startHandPos(1, 1) == modelParameters.initial(1, modelParameters.angle_n)
          noise(1, 1) = 0;
          scale(1, 1) = abs(test_data.startHandPos(1, 1) - modelParameters.initial(1, modelParameters.angle_n));
      elseif test_data.startHandPos(1, 1) > modelParameters.initial(1, modelParameters.angle_n)
          noise(1, 1) = 1;
          scale(1, 1) = abs(test_data.startHandPos(1, 1) - modelParameters.initial(1, modelParameters.angle_n));
      end
      
      if test_data.startHandPos(2, 1) < modelParameters.initial(2, modelParameters.angle_n)
          noise(2, 1) = -1;
          scale(2, 1) = abs(test_data.startHandPos(2, 1) - modelParameters.initial(2, modelParameters.angle_n));
      elseif test_data.startHandPos(2, 1) == modelParameters.initial(2, modelParameters.angle_n)
          noise(2, 1) = 0;
          scale(2, 1) = abs(test_data.startHandPos(2, 1) - modelParameters.initial(2, modelParameters.angle_n));
      elseif test_data.startHandPos(2, 1) > modelParameters.initial(2, modelParameters.angle_n)
          noise(2, 1) = 1;
          scale(2, 1) = abs(test_data.startHandPos(2, 1) - modelParameters.initial(2, modelParameters.angle_n));
      end
      modelParameters.noise = noise;
      modelParameters.scale = scale;
  end
  
  if time_point <= size(modelParameters.traces, 3)
      if time_point == 5
          
          templates = modelParameters.templates4;
          spikes = test_data.spikes;

          n_n = size(spikes, 1);
          n_bin = 400/100;
          spikes = spikes(:,1:400);

          B = movsum(spikes, 100, 2, 'Endpoints','discard');
          B = reshape(B(:,1:100:end), n_n*n_bin, 1);
          B = [mean(spikes, 2) ; B];
          B = [B(1:n_n) ./ repmat(sum( ...
              B(1:n_n)), n_n, 1); B];

          [~, angle] = min(sum(abs(templates - ...
                       repmat(B, 1, 8)' ...
                       ), 2));
          modelParameters.angle_n = mode(angle);
      end
      if time_point == 10
          
          templates = modelParameters.templates6;
          spikes = test_data.spikes;

          n_n = size(spikes, 1);
          n_bin = 500/125;
          spikes = spikes(:,1:500);

          B = movsum(spikes, 125, 2, 'Endpoints','discard');
          B = reshape(B(:,1:125:end), n_n*n_bin, 1);
          B = [mean(spikes, 2) ; B];
          B = [B(1:n_n) ./ repmat(sum( ...
              B(1:n_n)), n_n, 1); B];

          [~, angle] = min(sum(abs(templates - ...
                       repmat(B, 1, 8)' ...
                       ), 2));
          modelParameters.angle_n = mode(angle);
      end
      
      x = modelParameters.traces(modelParameters.angle_n, 1, time_point) + modelParameters.noise(1, 1)*modelParameters.deviation(modelParameters.angle_n, 1, time_point)*modelParameters.scale(1, 1)/7.5;
      y = modelParameters.traces(modelParameters.angle_n, 2, time_point) + modelParameters.noise(2, 1)*modelParameters.deviation(modelParameters.angle_n, 2, time_point)*modelParameters.scale(2, 1)/7.5;
  else
      x = modelParameters.objectives(1, modelParameters.angle_n);
      y = modelParameters.objectives(2, modelParameters.angle_n);
  end

end