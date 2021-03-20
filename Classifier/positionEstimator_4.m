function [angle, modelParameters] = positionEstimator_4(test_data, modelParameters)
  
  classifier = modelParameters.classifier;
  if size(test_data.spikes, 2) == 320
      angle = classifier.binSimilarityDistributions(modelParameters.templates, ...
                                                    modelParameters.distributions, ...
                                                    modelParameters.pref_neuron, ...
                                                    modelParameters.sum_activity, ...
                                                    test_data.spikes, 1, 300);
      modelParameters.angle_n = angle;
  end
  angle = modelParameters.angle_n;
 
end