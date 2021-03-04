%% Pre-Processing:

%% loading data
clc
clear all
load monkeydata_training.mat

%% classifying motor and pre-motor neurons
RMS_difference = zeros(size(trial(1, 1).spikes, 1), size(trial,1), size(trial,2));
for neuron=1:size(trial(1, 1).spikes,1)
    for angle_n=1:size(trial,2)
        for trial_n=1:size(trial,1)
            % length_ = size(trial(trial_n, angle_n).spikes,2);
            length_= 571;
            diff = mean(trial(trial_n, angle_n).spikes(neuron, 1:320)) - mean(trial(trial_n, angle_n).spikes(neuron, 321:length_));
            RMS_difference(neuron, trial_n, angle_n) = diff; 
        end
    end
end

neuron_n = 31;
angle_n = 1;
trial_n = 1;
plot(RMS_difference(neuron_n, :, 6))
yline(0)

which_neurons = zeros(size(trial(1, 1).spikes, 1), size(trial,2));
for neuron_n=1:size(trial(1, 1).spikes,1)
    for angle_n=1:size(trial,2)
        if length(find(RMS_difference(neuron_n, :, angle_n)<0)) < 20
            which_neurons(neuron_n, angle_n) = 1;
        end
    end
end

%% pre-processing for KNN:
neurons_selected = 1:98;
length_ = 300;
array_1d_neurons = zeros(size(trial,1)*size(trial,2), length(neurons_selected));
block_counter = 0;
for trial_n=1:size(trial,1)
    for angle_n=1:size(trial,2)
        block_counter = block_counter + 1;
        temp = zeros(1,length(neurons_selected));
        for neuron=1:length(neurons_selected)
            temp(neuron) = sum(trial(trial_n, angle_n).spikes(neuron, 1:length_));
        end
        array_1d_neurons(block_counter, :) = temp;
    end
end

hold on;
trials_ = 2:8:size(array_1d_neurons,1);
for i=1:length(trials_)
    plot(array_1d_neurons(trials_(1), :) - array_1d_neurons(trials_(i), :))
end

