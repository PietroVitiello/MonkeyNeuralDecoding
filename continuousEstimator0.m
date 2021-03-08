%%% Team Members: WRITE YOUR TEAM MEMBERS' NAMES HERE
%%% BMI Spring 2015 (Update 17th March 2015)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         PLEASE READ BELOW            %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function positionEstimator has to return the x and y coordinates of the
% monkey's hand position for each trial using only data up to that moment
% in time.
% You are free to use the whole trials for training the classifier.

% To evaluate performance we require from you two functions:

% A training function named "positionEstimatorTraining" which takes as
% input the entire (not subsampled) training data set and which returns a
% structure containing the parameters for the positionEstimator function:
% function modelParameters = positionEstimatorTraining(training_data)
% A predictor named "positionEstimator" which takes as input the data
% starting at 1ms and UP TO the timepoint at which you are asked to
% decode the hand position and the model parameters given by your training
% function:

% function [x y] = postitionEstimator(test_data, modelParameters)
% This function will be called iteratively starting with the neuronal data 
% going from 1 to 320 ms, then up to 340ms, 360ms, etc. until 100ms before 
% the end of trial.


% Place the positionEstimator.m and positionEstimatorTraining.m into a
% folder that is named with your official team name.

% Make sure that the output contains only the x and y coordinates of the
% monkey's hand.


function [modelParameters] = positionEstimatorTraining(training_data)
    % Arguments:

    % - training_data:
    %     training_data(n,k)              (n = trial id,  k = reaching angle)
    %     training_data(n,k).trialId      unique number of the trial
    %     training_data(n,k).spikes(i,t)  (i = neuron id, t = time)
    %     training_data(n,k).handPos(d,t) (d = dimension [1-3], t = time)

    % ... train your model

    % Return Value:

    % - modelParameters:
    %     single structure containing all the learned parameters of your
    %     model and which can be used by the "positionEstimator" function.
    
    %%%clean data for classifier
    silent_neuron = [8 10 11 38 49 52 73 74 76];
    trial = training_data;
    for angle_n = 1:size(training_data, 2)
        for trial_n = 1:size(training_data, 1)
            for neuron_n = 1:length(silent_neurons)
                trial(trial_n, angle_n).spikes(silent_neurons(neuron_n), :) = [];
                %trial(trial_n, angle_n).handPos(silent_neurons(neuron_n), :) = [];
            end
        end
    end
    
    %%%active neurons
    n = 4;
    pre_motor_window = 1 : 320;
    average_spike_trains = zeros(size(trial(1,1).spikes, 1), size(trial, 2));

    for angle_n = 1:size(trial, 2)
        for i = 1:size(trial, 1)
            average_spike_trains(:,angle_n) = average_spike_trains(:,angle_n) + mean(trial(i, angle_n).spikes(:, 1:pre_motor_window), 2);
        end
    end

    %active neurons is a matrix, each column represents one angle and
    %the neurons are ordered from the highest to lowest
    [~, all_active_neurons] = sort(average_spike_trains, 'descend');

    active_neurons = [];
    col = 1;
    row = 1;

    while length(active_neurons) < n*size(all_active_neurons, 2)
        temp = all_active_neurons(row, col);

        if isempty(find(active_neurons == temp, 1))
            active_neurons = [active_neurons temp];

            if length(active_neurons) == col*n
                col = col + 1;
                row = 0;
            end
        end

        row = row + 1;
    end
    
    
    %%% create matrix
    trial_struct = trial;
    training_percentage = 100;
    n_ms = 320;
    
    n_trials = size(trial_struct, 1);
    n_angles = size(trial_struct, 2);
    n_neurons = length(active_neurons);

    mx = zeros(n_trials, n_angles, n_neurons, n_ms);
    for trial = 1:n_trials
        for angle = 1:n_angles
            mx(trial,angle,:,:) = trial_struct(trial, angle).spikes(active_neurons, 1:n_ms);
        end
    end

    n_train_trials = round(n_trials*training_percentage / 100);
    train_matrix = zeros(n_train_trials, n_angles, n_neurons, n_ms);
    %test_mx = zeros(n_trials-n_train_trials, n_angles, n_neurons, n_ms);
    for angle = 1:n_angles
        rand_trials = randperm(n_trials)';
        train_matrix(:,angle,:,:) = mx(rand_trials(1:n_train_trials), angle, :, :);
        %test_mx(:,angle,:,:) = mx(rand_trials(n_train_trials+1:end), angle, :, :);
    end
            
          
    
    %%% MLE fitting
    [~, n_angles, n_neurons, ~] = size(train_matrix);

    train_data = sum(train_matrix, 4);
    par = zeros(n_neurons, n_angles, 2);

    for neuron = 1:n_neurons
        for angle = 1:n_angles
            par(neuron, angle, :) = mle(train_data(:, angle, neuron));
        end
    end
    
    modelParameters.classifier = par;
    
    %%% Trajectory average
    average_trj = zeros(8, 571);
    
    for a = 1:8
        for i = 1:size(trial, 1)
            for i = 1:571
                average_trj(a, :) = average_trj(a, :) + trial(






end

function [x, y] = positionEstimator(test_data, modelParameters)

  % **********************************************************
  %
  % You can also use the following function header to keep your state
  % from the last iteration
  %
  % function [x, y, newModelParameters] = positionEstimator(test_data, modelParameters)
  %                 ^^^^^^^^^^^^^^^^^^
  % Please note that this is optional. You can still use the old function
  % declaration without returning new model parameters. 
  %
  % *********************************************************

  % - test_data:
  %     test_data(m).trialID
  %         unique trial ID
  %     test_data(m).startHandPos
  %         2x1 vector giving the [x y] position of the hand at the start
  %         of the trial
  %     test_data(m).decodedHandPos
  %         [2xN] vector giving the hand position estimated by your
  %         algorithm during the previous iterations. In this case, N is 
  %         the number of times your function has been called previously on
  %         the same data sequence.
  %     test_data(m).spikes(i,t) (m = trial id, i = neuron id, t = time)
  %     in this case, t goes from 1 to the current time in steps of 20
  %     Example:
  %         Iteration 1 (t = 320):
  %             test_data.trialID = 1;
  %             test_data.startHandPos = [0; 0]
  %             test_data.decodedHandPos = []
  %             test_data.spikes = 98x320 matrix of spiking activity
  %         Iteration 2 (t = 340):
  %             test_data.trialID = 1;
  %             test_data.startHandPos = [0; 0]
  %             test_data.decodedHandPos = [2.3; 1.5]
  %             test_data.spikes = 98x340 matrix of spiking activity
  
  
  
  % ... compute position at the given timestep.
  
  % Return Value:
  
  % - [x, y]:
  %     current position of the hand
   
end

