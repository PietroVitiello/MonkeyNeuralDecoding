function [angle, modelParameters] = positionEstimator_9(test_data, modelParameters)
  
  classifier = modelParameters.classifier;
  if size(test_data.spikes, 2) == 320
      angle = classifier.findSimilarAngle_3D(modelParameters.templates1, test_data.spikes, 300, 1, 150);
      modelParameters.angle_n = angle;
  end
  angle = modelParameters.angle_n;
 
end