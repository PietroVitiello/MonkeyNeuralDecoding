% Test Script to give to the students, March 2015
%% Continuous Position Estimator Test Script
% This function first calls the function "positionEstimatorTraining" to get
% the relevant modelParameters, and then calls the function
% "positionEstimator" to decode the trajectory. 

function n_different_tests = testFunction_for_Classifier_6(teamName)

load monkeydata0.mat

n_different_tests = 20;
correct_first = 0;
correct_final = 0;
count_guess_error = zeros(1, 8);
count_correct_error = zeros(1, 8);
for test = 1:n_different_tests
    % Set random number generator
    %rng(2013);
    ix = randperm(length(trial));

    % Select training and testing data (you can choose to split your data in a different way if you wish)
    trainingData = trial(ix(1:80),:);
    testData = trial(ix(81:end),:);

    fprintf('Testing the classifier estimator...')
    
    n_predictions = size(testData,1)*8;

    % Train Model
    modelParameters = positionEstimatorTraining_6(trainingData);

    for tr=1:size(testData,1)
        display(['Decoding block ',num2str(tr),' out of ',num2str(size(testData,1))]);
        pause(0.001)
        for direc=randperm(8)

            times=320:20:size(testData(tr,direc).spikes,2);

            for t=times
                past_current_trial.trialId = testData(tr,direc).trialId;
                past_current_trial.spikes = testData(tr,direc).spikes(:,1:t);
                past_current_trial.startHandPos = testData(tr,direc).handPos(1:2,1); 

                if nargout('positionEstimator') == 2
                    [angle, newParameters] = positionEstimator_6(past_current_trial, modelParameters);
                    modelParameters = newParameters;
                elseif nargout('positionEstimator') == 1
                    [angle] = positionEstimator_6(past_current_trial, modelParameters);
                end

                if t == 320
                    correct_first = correct_first + (angle == direc);
                elseif t == times(end)
                    correct_final = correct_final + (angle == direc);
                end
                if angle ~= direc && t == 320
                    disp('�������������������')
                    for i = 1:8
                        if angle == i
                            count_guess_error(1, i) = count_guess_error(1, i) + 1;
                        end
                        if direc == i 
                            count_correct_error(1, i) = count_correct_error(1, i) + 1;
                        end
                    end
                end
            end
        end
    end
end
count_guess_error
count_correct_error

fprintf('\n\nCorrect predictions at the beginning: %.3f out of %d', correct_first/n_different_tests, n_predictions);

fprintf('\nCorrect final predictions: %.3f out of %d', correct_final/n_different_tests, n_predictions);

accuracy = correct_final / (n_predictions * n_different_tests);
fprintf('\nFinal accuracy: %f\n\n', accuracy);

rmpath(genpath(teamName))

end
