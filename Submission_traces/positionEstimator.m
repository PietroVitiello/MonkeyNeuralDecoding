function [x, y] = positionEstimator(test_data, modelParameters, direc)
  % Return Value:
  % - [x, y]:
  %     current position of the hand
  
  angle_n = direc;
  len = size(test_data.spikes,2);
  time_point = (len-300)/20;
  if time_point <= size(modelParameters.traces, 3)
      x = modelParameters.traces(angle_n, 1, time_point);
      y = modelParameters.traces(angle_n, 2, time_point);
  else
      x = modelParameters.objectives(1, angle_n);
      y = modelParameters.objectives(2, angle_n);
  end
 
end