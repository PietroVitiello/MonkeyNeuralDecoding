function [angle, modelParameters] = positionEstimator_2(test_data, modelParameters)
  
  classifier = modelParameters.classifier;
  if size(test_data.spikes, 2) == 320
      angle = classifier.findSimilarAngle2(modelParameters.templates1, test_data.spikes, 300, 1, 0.2);
%       angle = classifier.findSimilarAngle1_2(modelParameters.templates1, test_data.spikes, 300, 1);
      modelParameters.angle_n = angle;
%   elseif size(test_data.spikes, 2) == 400
%       angle = classifier.findSimilarAngle(modelParameters.templates2, test_data.spikes, 400, 1);
%       modelParameters.angle_n = angle;
  end
  angle = modelParameters.angle_n;
 
end