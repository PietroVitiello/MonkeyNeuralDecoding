function [angle, modelParameters] = positionEstimator_6(test_data, modelParameters)
 
  spikes = test_data.spikes;
  length_ = 320;
  classifier = modelParameters.classifier;
  angle_array = zeros(5, 1);
  if size(spikes, 2) <= length_
      init_spikes = sum(spikes(modelParameters.eccoli, 1:length_), 2);
      angle_array(1, 1) = predict(modelParameters.classifier1, init_spikes');
      
      angle_array(2, 1) = classifier.findSimilarAngle2(modelParameters.templates, test_data.spikes, 320, 1);
      
      angle_array(3, 1) = classifier.similarityDistributions(...
              modelParameters.templates, modelParameters.distributions ...
              , test_data.spikes, 1, 300);
      
      angle_array(4, 1) = classifier.apply_nDistribution_mle(test_data.spikes, ...
                                                 modelParameters.par, 1, 320);
                                
      init_spikes = sum(spikes(modelParameters.eccoli, 1:length_), 2);
      init_fr_diff =  zeros(1, size(modelParameters.templates, 1));
      for i = 1:size(modelParameters.templates, 1)
            init_fr_diff(1, i) = mean(abs(modelParameters.templates(i, :) - (sum(spikes(modelParameters.eccoli, 1:length_), 2)./length_)'));
      end
      init_spikes = [init_spikes' init_fr_diff.*20];
      angle_array(5, 1) = predict(modelParameters.classifier4, init_spikes);
      
      [most_freq num_occ] = mode(angle_array);

      modelParameters.angle_n = most_freq;
%       if (num_occ <= 2) && (size(unique(angle), 1) == 2) && (modelParameters.angle_n ~= angle(3, 1))
%           disp(modelParameters.angle_n);
%           disp(angle);
%           modelParameters.angle_n = angle(3, 1);
%           disp(modelParameters.angle_n);
%       end

  end
  
  angle = modelParameters.angle_n;
 
end