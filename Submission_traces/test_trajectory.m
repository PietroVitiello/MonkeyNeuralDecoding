clc
clear all
load monkeydata0
trj = Trajectory();
classifier = AngleClassifier();
process = Processing();
[trials, pos] = process.get_data_matrix(trial);

%% Test average length
lens = trj.movementDuration(trial);

%% Test average trajectory

templates = classifier.firingTemplate(trials);