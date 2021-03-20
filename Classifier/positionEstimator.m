function [angle, modelParameters] = positionEstimator(test_data, modelParameters)
 
  spikes = test_data.spikes;
  length_all = size(test_data.spikes,2);
  length_ = 320;
  classifier = modelParameters.classifier;
  if size(spikes, 2) <= length_
      init_spikes = sum(spikes(modelParameters.neurons, 1:length_), 2);
      %angle_n = predict(modelParameters.classifier1, init_spikes');
      angle_n = classifier.closestFiring(modelParameters.neuron_matrix, spikes, 1, 320);
      modelParameters.angle_n = angle_n;
  else
      if size(spikes, 2) == 360
          init_spikes = sum(spikes(modelParameters.eccoli, 1:360), 2);
          modelParameters.angle_n = predict(modelParameters.classifier2, init_spikes');
          %angle_n = classifier.closestFiring(modelParameters.neuron_matrix, spikes, 1, 360);
      end
      if size(spikes, 2) == 400
          init_spikes = sum(spikes(modelParameters.eccoli, 1:400), 2);
          modelParameters.angle_n = predict(modelParameters.classifier3, init_spikes');
      end
  end
  
  angle = modelParameters.angle_n;
 
end