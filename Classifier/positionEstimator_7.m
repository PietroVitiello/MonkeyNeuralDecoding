function [angle, modelParameters] = positionEstimator_7(test_data, modelParameters)
 
  spikes = test_data.spikes;
  length_ = 320;
  if size(spikes, 2) <= length_
      init_spikes = sum(spikes(modelParameters.neurons, 1:length_), 2);
      init_fr_diff =  zeros(1, size(modelParameters.templates1, 1));
      for i = 1:size(modelParameters.templates1, 1)
            init_fr_diff(1, i) = mean(abs(modelParameters.templates1(i, :) - (sum(spikes(modelParameters.neurons, 1:length_), 2)./length_)'));
      end
      init_spikes = [init_spikes' init_fr_diff.*20];
      angle = predict(modelParameters.knn, init_spikes);
      modelParameters.angle_n = angle;
  end
  angle = modelParameters.angle_n;
end