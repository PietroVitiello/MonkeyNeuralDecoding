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
trial_n = 100; %trial you want to observe
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

%% Raster per neuron across angles
min_length = 1*10^5;
for angle_n = 1:size(trial, 2)
    for i = 1:size(trial, 1)
        temp_length = size(trial(i,angle_n).spikes(1,:), 2);
        if temp_length < min_length
            min_length = temp_length;
        end
    end
end

% neurons = 1:14:size(trial(1, 1).spikes, 1);
neurons = 1:1:5;
overall_angles = zeros(length(neurons), min_length, size(trial, 2));
for angle_n=1:size(trial, 2)
    for neuron_n = 1:length(neurons)
        single_neuron = zeros(size(trial, 1), min_length);
        for trial_n=1:size(trial, 1)
            signal = trial(trial_n, angle_n).spikes(neurons(neuron_n), :);
            single_neuron(trial_n, :) = signal(1:min_length);
        end
        overall_angles(neuron_n, :, angle_n) = mean(single_neuron);
    end
end

max_ = -Inf;
for neuron_n=1:size(overall_angles,1)
    for angle_n=1:size(overall_angles,3)
        for time_point=1:size(overall_angles,2)
            if overall_angles(neuron_n, time_point, angle_n)> max_
                max_ = overall_angles(neuron_n, time_point, angle_n);
            end
        end
    end
end

c = 1;
for i=1:length(neurons)
    for j=1:size(trial, 2)
        figure(100)
        hold on;
        subplot(length(neurons),size(trial, 2), c)
        plot(1:min_length, smoothdata(smoothdata(overall_angles(i, :, j))));
        ylim([0 max_]);
        xline(300, '--r')
        c = c + 1;
    end
end

%% Linear fitting of labels
figure()
params_ = zeros(size(trial,2), 3);
for i=1:size(trial,2)
    temp = [];
    for j=1:size(trial,1)
        temp = [temp trial(j,i).handPos];
    end
    p = polyfit(temp(1, :), temp(2, :), 1);
    params_(i, :) = p;
    subplot(2,4,i)
    hold on;
    x = linspace(-20,100,1000);
    y = polyval(p, x);
    plot(x,y)
    scatter(temp(1, :), temp(2, :))
end

%% Hand Trajectory
angles = trial(1, :).handPos;
figure()
hold on
for i=1:size(trial,2)
    plot3(trial(1, i).handPos(1, :), trial(1, i).handPos(2, :), trial(1, i).handPos(3, :))
end
view(3)
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on

figure()
hold on
for i=1:size(trial,2)
    plot(trial(1, i).handPos(1, :), trial(1, i).handPos(2, :))
end
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on

%% Tuning curve
hold off;
pref_directions = zeros(1, size(trial(1,1).spikes, 1));
for unit_n=1:size(trial, 2)
    average_freq = zeros(1, size(trial, 2));
    std = zeros(1, size(trial, 2));
    angles = [30 70 110 150 190 230 310 350].*(pi/180); %angles corresponding to the 8

    for angle_n = 1:size(trial, 2)
        for i = 1:size(trial, 1)
            temp = average_freq(1, angle_n);
            x_new = mean(trial(i,angle_n).spikes(unit_n,:));
            average_freq(1, angle_n) = average_freq(1, angle_n) + x_new;
    %         std(1, angle_n) = 1/i*((i-1)*std(1, angle_n) + (i-1)*(i-2)*(average_freq(1, angle_n)-temp)^2); %running sample variance
            std(1, angle_n) = std(1, angle_n) + (x_new-temp)*(x_new-average_freq(1, angle_n)); %running sample intrmediate variance
        end
    end
    average_freq = average_freq/i; %averaging
    % std = sqrt(std); %standard deviation from variance
    std = sqrt(std/i); %standard deviation from intermediate variance
    pref_directions(unit_n) = find(average_freq == max(average_freq));
end 
%vertical bars representing the firing frequency at the different angles
figure()
hold on
bar(average_freq)
errorbar(1:8, average_freq, std/2, std/2)

%firing frequency in polar coordinates
figure()
polarplot([angles;angles], [zeros(size(average_freq));average_freq], 'LineWidth', 3);

%% averaged signal for each angle
min_length = 1*10^5;

for angle_n = 1:size(trial, 2)
    for i = 1:size(trial, 1)
        temp_length = size(trial(i,angle_n).spikes(1,:), 2);
        if temp_length < min_length
            min_length = temp_length;
        end
    end
end

average_spike_train = zeros(size(trial, 2), min_length);

for angle_n = 1:size(trial, 2)
    for i = 1:size(trial, 1)
        average_spike_train(angle_n, :) = average_spike_train(angle_n, :) + mean(trial(i,angle_n).spikes(:,1:min_length), 1);
    end
end
average_spike_train = average_spike_train/i;

for i = 1:size(trial, 2)
    figure(10)
    subplot(4,2,i)
    plot(average_spike_train(i,:));
    
    figure(11)
    subplot(4,2,i)
    %plot(smooth(smooth(smooth(average_spike_train(i,:)))));
end

%% Hand position delta
min_length = 1*10^5;
for angle_n = 1:size(trial, 2)
    for i = 1:size(trial, 1)
        temp_length = size(trial(i,angle_n).spikes(1,:), 2);
        if temp_length < min_length
            min_length = temp_length;
        end
    end
end

average_hand_pos = zeros(2, min_length);
average_deltaHand = zeros(size(trial, 2), min_length-1);

for angle_n = 1:size(trial, 2)
    for i = 1:size(trial, 1)
        average_hand_pos(1, :) = average_hand_pos(1, :) + trial(i,angle_n).handPos(1,1:min_length)/size(trial, 1);
        average_hand_pos(2, :) = average_hand_pos(2, :) + trial(i,angle_n).handPos(2,1:min_length)/size(trial, 1);
    end
%     average_deltaHand(angle_n, :) = pdist([average_hand_pos(angle_n, 2:end,:); average_hand_pos(angle_n, 1:end-1,:)], 'euclidean');
    for i = 1:min_length-1
        average_deltaHand(angle_n, i) = norm(average_hand_pos(:, i+1) - average_hand_pos(:, i));
    end
    average_hand_pos = zeros(2, min_length);
end

figure()
for i = 1:size(trial, 2)
    subplot(4,2,i)
    plot(average_deltaHand(i,:));
end

%% PSTH across angles
dummy = 10000;
for n_unit_i = 1:98
    for angle_i = 1:size(trial, 2)
        for trial_i = 1:size(trial, 1)
            n_bins = length(trial(trial_i, angle_i).spikes(n_unit_i, :));
            if n_bins < dummy
                dummy = n_bins;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%SUM ACROSS ANGLES%%%%%%%%%%%%%%%%
figure
count = 0;
n_neurons = 98;
step = 25;

useful_neurons = 1:98;

for n_unit_i = 17:32 %Choose neuron units to plot
    t = [];
    for angle_i = 1:size(trial, 2)
        for trial_i = 1:size(trial, 1)
            [~, t_i] = find(trial(trial_i, angle_i).spikes(n_unit_i, 1:dummy) > 0);
            t = [t t_i];
        end     
    end
    subplot(4, 4, count+1)
    hold on
    histogram(t, dummy)
    xline(300, 'color', 'b');
    xlim([0 dummy])
    ylim([0 100])
    [F, xi] = ksdensity(t);
    plot(xi, F*length(t), 'LineWidth', 2)
    smooth = F*length(t);
    if (max(smooth) < 15)
        useful_neurons(1, n_unit_i) = 0;
    end
    
    count = count + 1;
end

useless_neurons = find(useful_neurons == 0);

%% Highest response per angle
pre_motor_window = 320;

average_spike_trains = zeros(size(trial(1,1).spikes, 1), size(trial, 2));

for angle_n = 1:size(trial, 2)
    for i = 1:size(trial, 1)
        average_spike_trains(:,angle_n) = average_spike_trains(:,angle_n) + mean(trial(i, angle_n).spikes(:, 1:pre_motor_window), 2);
    end
end

%active neurons is a matrix, each column represents one angle and
%the neurons are ordered from the highest to lowest
[~, active_neurons] = sort(average_spike_trains, 'descend');
