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


%% PSTH
trial_n = 1; %trial you want to observe
angle_n = 1; %reaching angle you want ot observe
bin_n = length(trial(trial_n,angle_n).spikes(1,:)); %number of bins equal
                                                    %to the # of ms

[~, col] = find(trial(trial_n,angle_n).spikes(:,:) > 0);
histogram(col, length(trial(trial_n,angle_n).spikes(1,:)))
normalization_factor = length(col); %factor to calculate the frequency from the probability
histogram(col, bin_n)

hold on
[f,xi] = ksdensity(col); %calculates probability distribution
plot(xi,f*normalization_factor, 'LineWidth', 2)

%% Hand Trajectory
angles = trial(1, :).handPos;
figure(2)
hold on
for i=1:size(trial,2)
    plot3(trial(1, i).handPos(1, :), trial(1, i).handPos(2, :), trial(1, i).handPos(3, :))
end
view(3)
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on

figure(3)
hold on
for i=1:size(trial,2)
    plot(trial(1, i).handPos(1, :), trial(1, i).handPos(2, :))
end
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on

%% Tuning curve
unit_n = 1;
average_freq = zeros(1, size(trial, 2));
angles = [30 70 110 150 190 230 310 350].*(pi/180);

for angle_n = 1:size(trial, 2)
    for i = 1:size(trial, 1)
        average_freq(1, angle_n) = average_freq(1, angle_n) + mean(trial(i,angle_n).spikes(unit_n,:));
    end
end
average_freq = average_freq/i;
bar(average_freq)

polarplot([angles;angles], [zeros(size(average_freq));average_freq], 'LineWidth', 3);
