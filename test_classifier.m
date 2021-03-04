clc
clear all
load monkeydata_training

processor = Processing();

active_neurons = processor.mostActive(trial, 4)