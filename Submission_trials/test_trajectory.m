clc
clear all
load monkeydata0
trj = Trajectory();
classifier = AngleClassifier();

%% Test average length
lens = trj.movementDuration(trial);

%% Test average trajectory

trajectories = classifier.meanTraces(trial, 1, 320);