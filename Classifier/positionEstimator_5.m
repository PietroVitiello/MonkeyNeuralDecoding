function [angle, modelParameters] = positionEstimator_5(test_data, modelParameters)
  
  classifier = modelParameters.classifier;
  if size(test_data.spikes, 2) == 320
      angle = classifier.apply_nDistribution_mle(test_data.spikes, ...
                                                 modelParameters.par, 1, 300);
      modelParameters.angle_n = angle;
  end
  angle = modelParameters.angle_n;
 
end