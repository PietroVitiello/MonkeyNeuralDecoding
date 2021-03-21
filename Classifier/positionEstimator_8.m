function [angle, modelParameters] = positionEstimator_8(test_data, modelParameters, direc)
  
  classifier = modelParameters.classifier;
  if size(test_data.spikes, 2) == 320
      angle = classifier.SimilarityMLE(modelParameters.templates1, modelParameters.par, test_data.spikes, 300, 1);
      modelParameters.angle_n = angle;
%       [vector, angle] = classifier.visualizeSimilarityVectors(modelParameters.templates1, test_data.spikes, 300, 1);
%       modelParameters.angle_n = angle;
%       if angle == direc
%           vector
%       else
%           vector'
%       end
%       vector(direc)
  end
  angle = modelParameters.angle_n;
 
end