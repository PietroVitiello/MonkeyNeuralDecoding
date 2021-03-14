% Test Script to give to the students, March 2015
%% Continuous Position Estimator Test Script
% This function first calls the function "positionEstimatorTraining" to get
% the relevant modelParameters, and then calls the function
% "positionEstimator" to decode the trajectory. 

function RMSE = testFunction_for_students_MTb(teamName)

load monkeydata0.mat

% Set random number generator
rng(2013);
ix = randperm(length(trial));

addpath(teamName);

% Select training and testing data (you can choose to split your data in a different way if you wish)
trainingData = trial(ix(1:50),:);
testData = trial(ix(51:end),:);

fprintf('Testing the classifier estimator...')

correct_first = 0;
correct_final = 0;
n_predictions = size(testData,1)*8;

% Train Model
modelParameters = positionEstimatorTraining(trainingData);

for tr=1:size(testData,1)
    display(['Decoding block ',num2str(tr),' out of ',num2str(size(testData,1))]);
    pause(0.001)
    for direc=randperm(8)

        times=320:20:size(testData(tr,direc).spikes,2);
        
        for t=times
            past_current_trial.trialId = testData(tr,direc).trialId;
            past_current_trial.spikes = testData(tr,direc).spikes(:,1:t);

            past_current_trial.startHandPos = testData(tr,direc).handPos(1:2,1); 
            
            if nargout('positionEstimator') == 3
                [angle, newParameters] = positionEstimator(past_current_trial, modelParameters);
                modelParameters = newParameters;
            elseif nargout('positionEstimator') == 2
                [angle] = positionEstimator(past_current_trial, modelParameters);
            end
            
            if t == 320
                correct_first = correct_first + (angle == direc);
            elseif t == times(end)
                correct_final = correct_final + (angle == direc);
            end
        end
    end
end

sprintf('Correct predictions from first 320: %d out of %d', correct_first, n_predictions);

sprintf('Correct final predictions: %d out of %d', correct_final, n_predictions);

accuracy = correct_final / n_predictions;
sprintf('Final accuracy: %f', accuracy)

rmpath(genpath(teamName))

end
