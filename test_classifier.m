clc
clear all
load monkeydata_training

silent_neuron = [8 10 11 38 49 52 73 74 76];
processor = Processing();
clean_trial = processor.clean_dataset(trial, silent_neuron);

active_neurons = processor.mostActive(clean_trial, 4);