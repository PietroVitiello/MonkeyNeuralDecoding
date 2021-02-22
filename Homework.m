clc
clear all
load monkeydata_training.mat

%% Raster plot same trial different unit
[row, col] = find(trial(1,1).spikes(:,:) > 0);
scatter(col, row, 2, 'filled')
%spikes(:,:) = find(trial(1,1).spikes(:,:) > 0);
%scatter(spikes, ones(size(spikes))*[1:size(spikes,2)])

%plotSpikeRaster(trial(1,1).spikes(1,:))

%% Raster plot same unit different trials
row = [];
col = [];
for i = 1:size(trial, 1)
    [~, col_temp] = find(trial(i,1).spikes(1,:) > 0);
    row = [row ones(1,length(col_temp))*i];
    col = [col col_temp];
end
%row = ones(size(pos))*[1:size(pos, 1)]
scatter(col, row, 2, 'filled')
