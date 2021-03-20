function [angle, modelParameters] = positionEstimator_3(test_data, modelParameters)
  
  classifier = modelParameters.classifier;
  if size(test_data.spikes, 2) == 320
      angle = classifier.similarityDistributions(...
              modelParameters.templates, modelParameters.distributions ...
              , test_data.spikes, 1, 300);
      modelParameters.angle_n = angle;
  end
  angle = modelParameters.angle_n;
 
end