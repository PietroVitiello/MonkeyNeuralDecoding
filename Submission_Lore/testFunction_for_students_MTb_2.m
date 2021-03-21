% Test Script to give to the students, March 2015
%% Continuous Position Estimator Test Script
% This function first calls the function "positionEstimatorTraining" to get
% the relevant modelParameters, and then calls the function
% "positionEstimator" to decode the trajectory. 

function [average_RMSE, n_different_tests] = testFunction_for_students_MTb_2(teamName)

load monkeydata0.mat

n_different_tests = 20;
RMSE_runs = zeros(1, n_different_tests);
correct_first = 0;
correct_final = 0;

for test = 1:n_different_tests
    % Set random number generator
    %rng(2013);
    ix = randperm(length(trial));

    % Select training and testing data (you can choose to split your data in a different way if you wish)
    trainingData = trial(ix(1:80),:);
    testData = trial(ix(81:end),:);

    fprintf('Testing the classifier estimator...')
    
    meanSqError = 0;
    n_predictions = 0;
    
    figure
    hold on
    axis square
    grid

    % Train Model
    modelParameters = positionEstimatorTraining_2(trainingData);

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
                    [decodedPosX, decodedPosY, newParameters] = positionEstimator_2(past_current_trial, modelParameters);
                    modelParameters = newParameters;
                elseif nargout('positionEstimator') == 2
                    [decodedPosX, decodedPosY] = positionEstimator_2(past_current_trial, modelParameters);
                end
                
                decodedPos = [decodedPosX; decodedPosY];
                decodedHandPos = [decodedHandPos decodedPos];
                meanSqError = meanSqError + norm(testData(tr,direc).handPos(1:2,t) - decodedPos)^2;

%                 if t == 320
%                     correct_first = correct_first + (angle == direc);
%                 elseif t == times(end)
%                     correct_final = correct_final + (angle == direc);
%                 end
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

% fprintf('\n\nCorrect predictions at the beginning: %.3f out of %d', correct_first/n_different_tests, n_predictions);
% 
% fprintf('\nCorrect final predictions: %.3f out of %d', correct_final/n_different_tests, n_predictions);
% 
% accuracy = correct_final / (n_predictions * n_different_tests);
% fprintf('\nFinal accuracy: %f\n\n', accuracy);
% 
% rmpath(genpath(teamName))

end