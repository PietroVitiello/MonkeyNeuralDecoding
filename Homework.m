clc
clear all
load monkeydata_training.mat

%% Raster plot

[row, col] = find(trial(1,1).spikes(:,:) > 0);
scatter(row, col, 2, 'filled')
%spikes(:,:) = find(trial(1,1).spikes(:,:) > 0);
%scatter(spikes, ones(size(spikes))*[1:size(spikes,2)])

%plotSpikeRaster(trial(1,1).spikes(1,:))