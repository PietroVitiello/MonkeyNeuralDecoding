function [x, y, modelParameters] = positionEstimator_4(test_data, modelParameters)
  % Return Value:
  % - [x, y]:
  %     current position of the hand
  
  classifier = modelParameters.classifier;
  
  len = size(test_data.spikes,2);
  time_point = (len-300)/20;
  noise = zeros(2, 1);
  scale = zeros(2, 1);
  angle = zeros(1, 1);
  
  if time_point == 1
%       angle(1, 1) = classifier.similarityDistributions(...
%               modelParameters.templates, modelParameters.distributions ...
%               , test_data.spikes, 1, 300);
      
      %init_spikes = sum(test_data.spikes(1:98, 1:len), 2);
      %angle(1, 1) = predict(modelParameters.classifier1, init_spikes');
      
      %angle(2, 1) = classifier.findSimilarAngle2(modelParameters.templates, test_data.spikes, 300, 1);
      
%       angle(2, 1) = classifier.findSimilarAngle_3D(modelParameters.templates1, test_data.spikes, 300, 1, 150);
      
      angle(1, 1) = classifier.findSimilarAngle_2n3D(modelParameters.templates2, test_data.spikes, 300, 1, 100);
      
      modelParameters.angle_n = mode(angle);
      
%       if (angle(1, 1) ~= angle(2, 1)) || (angle(1, 1) ~= angle(3, 1)) || (angle(2, 1) ~= angle(3, 1))
%           modelParameters.angle_n = angle(1, 1);
%           disp('################################"')
%       end
      
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
%       if time_point == 3
%           angle(1, 1) = classifier.findSimilarAngle_2n3D(modelParameters.templates3, test_data.spikes, len, 1, 120);
%           modelParameters.angle_n = mode(angle);
%       end
      if time_point == 5
          angle(1, 1) = classifier.findSimilarAngle_2n3D(modelParameters.templates4, test_data.spikes, len, 1, 100);
          modelParameters.angle_n = mode(angle);
      end
%       if time_point == 7
%           angle(1, 1) = classifier.findSimilarAngle_2n3D(modelParameters.templates5, test_data.spikes, len, 1, 110);
%           modelParameters.angle_n = mode(angle);
%       end
      if time_point == 10
          angle(1, 1) = classifier.findSimilarAngle_2n3D(modelParameters.templates6, test_data.spikes, len, 1, 125);
          modelParameters.angle_n = mode(angle);
      end
      
      x = modelParameters.traces(modelParameters.angle_n, 1, time_point) + modelParameters.noise(1, 1)*modelParameters.deviation(modelParameters.angle_n, 1, time_point)*modelParameters.scale(1, 1)/7.5;
      y = modelParameters.traces(modelParameters.angle_n, 2, time_point) + modelParameters.noise(2, 1)*modelParameters.deviation(modelParameters.angle_n, 2, time_point)*modelParameters.scale(2, 1)/7.5;
  else
      x = modelParameters.objectives(1, modelParameters.angle_n);
      y = modelParameters.objectives(2, modelParameters.angle_n);
  end

end