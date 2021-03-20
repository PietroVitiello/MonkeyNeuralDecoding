% Test Script to give to the students, March 2015
%% Continuous Position Estimator Test Script
% This function first calls the function "positionEstimatorTraining" to get
% the relevant modelParameters, and then calls the function
% "positionEstimator" to decode the trajectory. 

    function testFunction_for_students_MTb(teamName)

    load monkeydata0.mat
    
    n_different_tests = 10;
    training_time = 0;
    testing_time = 0;
    meanSqError = 0;
    n_predictions = 0;
    
    figure
    hold on
    axis square
    grid

    for test = 1:n_different_tests
        % Select training and testing data (you can choose to split your data in a different way if you wish)
        % Set random number generator
        % rng(2013);
        ix = randperm(length(trial));
        trainingData = trial(ix(1:80),:);
        testData = trial(ix(81:end),:);

        fprintf('Testing the continuous position estimator...')

        % Train Model
        tic
        modelParameters = positionEstimatorTraining(trainingData);
        training_time = training_time + toc;

        tic
        for tr=1:size(testData,1)
            display(['Decoding block ',num2str(tr),' out of ',num2str(size(testData,1))]);
            pause(0.001)
            for direc=randperm(8) 
                decodedHandPos = [];

                times=320:20:size(testData(tr,direc).spikes,2);
        %         times = 320;
                for t=times
                    past_current_trial.trialId = testData(tr,direc).trialId;
                    past_current_trial.spikes = testData(tr,direc).spikes(:,1:t); 
                    past_current_trial.decodedHandPos = decodedHandPos;

                    past_current_trial.startHandPos = testData(tr,direc).handPos(1:2,1); 

                    if nargout('positionEstimator') == 3
                        [decodedPosX, decodedPosY, newParameters] = positionEstimator(past_current_trial, modelParameters);
                        modelParameters = newParameters;
                    elseif nargout('positionEstimator') == 2
                        [decodedPosX, decodedPosY] = positionEstimator(past_current_trial, modelParameters);
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
        testing_time = testing_time + toc;
    end

    legend('Decoded Position', 'Actual Position')

    RMSE = sqrt(meanSqError/(n_predictions))

    avg_train_time = training_time/n_different_tests;
    avg_test_time = testing_time/n_different_tests;
    fprintf('\nThe average elapsed time for training is %f', avg_train_time);
    fprintf('\nThe average elapsed time for testing is %f', avg_test_time);
    fprintf('\nThe estimated total time for 80 trials would have been %f\n\n', ...
            avg_train_time + avg_test_time*4);

    end
