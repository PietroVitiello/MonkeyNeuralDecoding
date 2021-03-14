clc
load monkeydata0
testFunction_for_students_MTb('Team Cassius')


% estimator = PositionEstimator_cl();
% processor = Processing();
% 
% neurons_per_angle = 5;
% mode = 'vector';
% if strcmp(mode, 'matrix')
%   neurons_estimator_matrix = processor.mostActive(trial, neurons_per_angle, 300, 571, mode);
% end
% if strcmp(mode, 'vector')
%   neurons_estimator = processor.mostActive(trial, neurons_per_angle, 300, 571, mode); 
%   neurons_estimator_matrix = zeros(length(neurons_estimator), size(trial,2));
%   for angle_n = 1:size(trial, 2)
%       neurons_estimator_matrix(:, angle_n) = neurons_estimator';
%   end
% end
% 
% [state0, eeg_train, ~, x_train, ~] = estimator.getDataset(trial, 1, 100);
% eeg_train = estimator.non_redundant(eeg_train, neurons_estimator_matrix);
% [eeg_train_, x_train_] = processor.get_average_data(eeg_train, x_train);
