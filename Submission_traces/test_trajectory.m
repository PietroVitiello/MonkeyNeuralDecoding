clc
clear all
load monkeydata0
trj = Trajectory();
classifier = AngleClassifier();
process = Processing();

%% Test average length
lens = trj.movementDuration(trial);

% %% Test average trajectory
% trajectories = classifier.meanTraces(trial, 1, 320);

%% Create matrices
[trials, pos] = process.get_data_matrix(trial);

%% Find objectives
obj = trj.objective_positions(pos);

figure()
scatter(obj(1,:), obj(2,:))

%% Average of trajectories
[avgT, stdT] = trj.averageTrajectory(pos);

% figure()
hold on
for i = 1:8
    scatter(avgT(i,1,:), avgT(i,2,:))
end


