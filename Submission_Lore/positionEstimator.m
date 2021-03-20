function [x, y, modelParameters] = positionEstimator(test_data, modelParameters)
  % Return Value:
  % - [x, y]:
  %     current position of the hand
  
  classifier = modelParameters.classifier;
  
  len = size(test_data.spikes,2);
  time_point = (len-300)/20;
  noise = zeros(2, 1);
  scale = zeros(2, 1);
  
  if time_point == 1
      angle = classifier.similarityDistributions(...
              modelParameters.templates, modelParameters.distributions ...
              , test_data.spikes, 1, 300);
      modelParameters.angle_n = angle;
      
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
      elseif test_data.startHandPos(2, 1) == modelParameters.initial(2, angle_n)
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
      x = modelParameters.traces(modelParameters.angle_n, 1, time_point) + modelParameters.noise(1, 1)*modelParameters.deviation(modelParameters.angle_n, 1, time_point)*modelParameters.scale(1, 1)/7.5;
      y = modelParameters.traces(modelParameters.angle_n, 2, time_point) + modelParameters.noise(2, 1)*modelParameters.deviation(modelParameters.angle_n, 2, time_point)*modelParameters.scale(2, 1)/7.5;
  else
      x = modelParameters.objectives(1, modelParameters.angle_n);
      y = modelParameters.objectives(2, modelParameters.angle_n);
  end
 
end