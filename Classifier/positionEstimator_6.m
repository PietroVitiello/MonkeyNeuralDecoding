function [angle, modelParameters] = positionEstimator_6(test_data, modelParameters)
 
  spikes = test_data.spikes;
  length_ = 320;
  classifier = modelParameters.classifier;
  angle = zeros(4, 1);
  if size(spikes, 2) <= length_
      init_spikes = sum(spikes(modelParameters.neurons, 1:length_), 2);
      angle(1, 1) = predict(modelParameters.classifier1, init_spikes');
      
      angle(2, 1) = classifier.findSimilarAngle2(modelParameters.templates, test_data.spikes, 300, 1);
      
      angle(3, 1) = classifier.similarityDistributions(...
              modelParameters.templates, modelParameters.distributions ...
              , test_data.spikes, 1, 300);
      
      angle(4, 1) = classifier.apply_nDistribution_mle(test_data.spikes, ...
                                                 modelParameters.par, 1, 300);
      
      [most_freq num_occ] = mode(angle);

      modelParameters.angle_n = most_freq;
      disp(angle);
%   else
%       if size(spikes, 2) == 360
%           init_spikes = sum(spikes(modelParameters.eccoli, 1:360), 2);
%           modelParameters.angle_n = predict(modelParameters.classifier2, init_spikes');
%           %angle_n = classifier.closestFiring(modelParameters.neuron_matrix, spikes, 1, 360);
%       end
%       if size(spikes, 2) == 400
%           init_spikes = sum(spikes(modelParameters.eccoli, 1:400), 2);
%           modelParameters.angle_n = predict(modelParameters.classifier3, init_spikes');
%       end
  end
  
  angle = modelParameters.angle_n;
 
end