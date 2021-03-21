% Test Script to give to the students, March 2015
%% Continuous Position Estimator Test Script
% This function first calls the function "positionEstimatorTraining" to get
% the relevant modelParameters, and then calls the function
% "positionEstimator" to decode the trajectory. 

function [average_RMSE, n_different_tests] = testFunction_for_students_MTb_3(teamName)

load monkeydata0.mat

n_different_tests = 20;
RMSE_runs = zeros(1, n_different_tests);

for test = 1:n_different_tests
    % Set random number generator
    %rng(2013);
    ix = randperm(length(trial));
    % addpath(teamName);

    % Select training and testing data (you can choose to split your data in a different way if you wish)
    trainingData = trial(ix(1:80),:);
    testData = trial(ix(81:end),:);

    fprintf('Testing the continuous position estimator...')

    meanSqError = 0;
    n_predictions = 0;

    figure
    hold on
    axis square
    grid

    % Train Model
    modelParameters = positionEstimatorTraining_3(trainingData);

    for tr=1:size(testData,1)
        display(['Decoding block ',num2str(tr),' out of ',num2str(size(testData,1))]);
        pause(0.001)
        for direc=randperm(8) 
            decodedHandPos = [];
            times=320:20:size(testData(tr,direc).spikes,2);
            for t=times
                past_current_trial.trialId = testData(tr,direc).trialId;
                past_current_trial.spikes = testData(tr,direc).spikes(:,1:t); 
                past_current_trial.decodedHandPos = decodedHandPos;

                past_current_trial.startHandPos = testData(tr,direc).handPos(1:2,1); 

                if nargout('positionEstimator') == 3
                    [decodedPosX, decodedPosY, newParameters] = positionEstimator_3(past_current_trial, modelParameters);
                    modelParameters = newParameters;
                elseif nargout('positionEstimator') == 2
                    [decodedPosX, decodedPosY] = positionEstimator_3(past_current_trial, modelParameters);
                end

                decodedPos = [decodedPosX; decodedPosY];
                decodedHandPos = [decodedHandPos decodedPos];

                meanSqError = meanSqError + norm(testData(tr,direc).handPos(1:2,t) - decodedPos)^2;

            end
            n_predictions = n_predictions+length(times);
            hold on
            plot(decodedHandPos(1,:),decodedHandPos(2,:), 'r');
            plot(testData(tr,direc).handPos(1,times),testData(tr,direc).handPos(2,times),'b')
        end
    end

    legend('Decoded Position', 'Actual Position')

    RMSE = sqrt(meanSqError/n_predictions);
    RMSE_vector(1, test) = RMSE;
    
end

average_RMSE = mean(RMSE_vector);
fprintf('\nAverage RMSE: %f\n\n', average_RMSE);

rmpath(genpath(teamName))

end