classdef Processing
    
    properties
      
    end
    
    methods
        function [trial, neurons] = clean_dataset(~, trial, silent_neurons)
            for angle_n = 1:size(trial, 2)
                for trial_n = 1:size(trial, 1)
                    trial(trial_n, angle_n).spikes(silent_neurons, :) = [];
                end
            end
            neurons = 1:98;
            neurons(silent_neurons) = [];
        end
        
        
        function [active_neurons, active_neuron_matrix] = mostActive(~, trial, n, neurons, lower_bound, upper_bound)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            The purpose of this function is to create a list of neurons
            that have a high activity in the first 320ms to use for the 
            angle classifier
            
            -input
            trial: is the data as a matrix of structs
            n: number of most active neurons taken for each angle without
               repetition
            lower_bound: lower bound of the sample window to consider
            upper_bound: upper bound of the sample window to consider
            
            -output
            active_neurons: 1 x n*(number of angles) list of neurons
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
            if nargin < 6
                lower_bound = 1;
            end
            if nargin < 5
                upper_bound = 320;
            end
            
            pre_motor_window = lower_bound : upper_bound;

            average_spike_trains = zeros(size(trial(1,1).spikes, 1), size(trial, 2));

            for angle_n = 1:size(trial, 2)
                for i = 1:size(trial, 1)
                    average_spike_trains(:,angle_n) = average_spike_trains(:,angle_n) + mean(trial(i, angle_n).spikes(:, 1:pre_motor_window), 2);
                end
            end

            %active neurons is a matrix, each column represents one angle and
            %the neurons are ordered from the highest to lowest
            [~, all_active_neurons] = sort(average_spike_trains, 'descend');
            all_active_neurons;
            
            active_neuron_matrix = zeros(n, size(trial, 2));
            for angle_n = 1:size(active_neuron_matrix, 2)
                active_neuron_matrix(:, angle_n) = neurons(all_active_neurons(1:n, angle_n));
            end
            
            active_neurons = [];
            col = 1;
            row = 1;
            
            while length(active_neurons) < n*size(all_active_neurons, 2)
                temp_index = all_active_neurons(row, col);
                temp = neurons(temp_index);
                
                if isempty(find(active_neurons == temp, 1))
                    active_neurons = [active_neurons temp];
                    
                    if col == size(all_active_neurons, 2)
                        row = row + 1;
                        col = 0;
                    end
                else
                    if col == size(all_active_neurons, 2)
                        row = row + 1;
                        col = 0;
                    end
                end
                
                col = col + 1;
            end 
        end
        
        function most_active = mostActiveV2(~, data_matrix)
            most_active = zeros(size(data_matrix,3), size(data_matrix,2));
            diffs = zeros(size(data_matrix, 3), 3);
            mean_fr = zeros(size(data_matrix,3), size(data_matrix,1));
            for neuron_n = 1:size(data_matrix,3)
                for trial_n = 1:size(data_matrix,1)
                    for angle_n = 1:size(data_matrix,2)
                        mean_fr(neuron_n, angle_n) = mean_fr(neuron_n, angle_n) + sum(data_matrix(trial_n, angle_n, neuron_n, 1:300));
                    end
                end
                mean_fr(neuron_n, :) = mean_fr(neuron_n, :)./size(data_matrix,1);
                [fr_sorted, idxs] = sort(mean_fr(neuron_n, :), 'descend');
                diffs(neuron_n, 1) = (fr_sorted(1) - fr_sorted(2))/fr_sorted(1);
                diffs(neuron_n, 2) = neuron_n; 
                diffs(neuron_n, 3) = idxs(1);
            end
            diffs = sortrows(diffs, 3);
            idx = 1;
            for angle_n = 1:size(data_matrix,2)
                count = 0;
                while diffs(idx,3) == angle_n 
                    count = count + 1;
                    if idx < size(diffs, 1)
                        idx = idx + 1;
                    else
                        break
                    end
                end
                temp = diffs(idx - count: idx-1, 1:2);
                temp = sortrows(temp, 1);
                most_active(1:size(temp,1), angle_n) = temp(:, 2);
            end
        end
        
        function [samples, labels] = create_dataset(~, trial, active_neurons, length_premotor, start_premotor)
            if nargin < 4
                start_premotor = 1;
            end
            if nargin < 3
                length_premotor = 320;
            end
            
            dataset = zeros(size(trial,1)*size(trial,2), length(active_neurons)+1);
            traj_count = 0;
            for angle_n = 1:size(trial,2)
                for trial_n = 1:size(trial, 1)
                    traj_count = traj_count + 1;
                    temp = zeros(1, length(active_neurons));
                    for neuron_n = 1:length(active_neurons)
                        temp(neuron_n) = sum(trial(trial_n, angle_n).spikes(active_neurons(neuron_n), start_premotor:length_premotor));
                        temp(length(active_neurons)+1) = angle_n;
                    end
                dataset(traj_count, :) = temp;
                end
            end
            % rand_d = dataset(randperm(size(dataset, 1)), :);
            samples = dataset(:, 1:length(active_neurons));
            labels = dataset(:, length(active_neurons)+1);
        end
        
        
        function [train_mx, test_mx] = data_as_matrix(~, trial_struct, active_neurons, training_percentage, n_ms)
            
            if nargin < 5
                n_ms = 320;
            end
            
            n_trials = size(trial_struct, 1);
            n_angles = size(trial_struct, 2);
            n_neurons = length(active_neurons);
            
            mx = zeros(n_trials, n_angles, n_neurons, n_ms);
            for trial = 1:n_trials
                for angle = 1:n_angles
                    mx(trial,angle,:,:) = trial_struct(trial, angle).spikes(active_neurons, 1:n_ms);
                end
            end
            
            n_train_trials = round(n_trials*training_percentage / 100);
            train_mx = zeros(n_train_trials, n_angles, n_neurons, n_ms);
            test_mx = zeros(n_trials-n_train_trials, n_angles, n_neurons, n_ms);
            for angle = 1:n_angles
                rand_trials = randperm(n_trials)';
                train_mx(:,angle,:,:) = mx(rand_trials(1:n_train_trials), angle, :, :);
                test_mx(:,angle,:,:) = mx(rand_trials(n_train_trials+1:end), angle, :, :);
            end
        end
        
        
        function mx = switchDim(~, input)
            
            [d1, d2, d3] = size(input);
            mx = zeros(d1, d3, d2);

            for i = 1:d1
                for j = 1:d3
                    for k = 1:d2
                        mx(i,j,k) = input(i,k,j);
                    end
                end
            end
            
        end
        
        
        function cov_mx = covariance(~, data_mx, angle)
            
            n_neurons = size(data_mx, 3);
            
            mx = sum(data_mx, 4);
            mx = reshape(mx(:,angle,:), [], n_neurons);
            
            cov_mx = cov(mx);
            
        end
        
        function [spikes_matrix, labels_matrix]  = get_data_matrix(~, trial)
            min_length = 571;
            spikes_matrix = zeros(size(trial,1), size(trial,2), size(trial(1,1).spikes, 1), min_length);
            labels_matrix = zeros(size(trial,1), size(trial,2), size(trial(1,1).handPos, 1)-1, min_length);
            for trial_n = 1:size(trial,1)
                for angle_n = 1: size(trial,2)
                    spikes_matrix(trial_n, angle_n, :, :) = trial(trial_n, angle_n).spikes(:, 1:min_length);
                    labels_matrix(trial_n, angle_n, :, :) = trial(trial_n, angle_n).handPos(1:2, 1:min_length);
                end
            end
        end
        
        function [firing_rates, angles] = get_dataset(~, templates, data_matrix, neurons_selected, start, stop)
            if nargin == 4
                start = 1;
                stop = 300;
            end
            
            if nargin < 4
                start = 1;
                stop = 300;
                neurons_selected = 1:size(data_matrix, 3);
            end
            
            firing_rates = zeros(size(data_matrix, 1)*size(data_matrix, 2)*size(data_matrix, 2), length(neurons_selected));
            angles = zeros(size(data_matrix, 1)*size(data_matrix, 2), 1);
            idx = 0;
            for angle_n = 1:size(data_matrix, 2)
                for trial_n = 1:size(data_matrix, 1)
                    for i = 1:1:size(data_matrix, 2)
                        idx = idx + 1;
                        firing_rates(idx, :) = templates(i, :) - (sum(squeeze(data_matrix(trial_n, angle_n, neurons_selected, start:stop)), 2)./(stop-start+1))';
                        angles(idx) = angle_n;
                    end
                end
            end
            temp = [firing_rates angles];
            rand_d = temp(randperm(size(temp, 1)), :);
            firing_rates = rand_d(:, 1:size(firing_rates, 2));
            angles = rand_d(:, size(firing_rates, 2)+1);
            
        end
        
        function av_diff = get_dataset_2(~, templates, data_matrix, neurons_selected, stop, start)
            if nargin < 5 
                start = 1;
                stop = 300;
            end
            
            if nargin < 4
                neurons_selected = 1:size(data_matrix, 3);
            end
            
            av_diff = zeros(size(data_matrix, 1)*size(data_matrix, 2), size(data_matrix, 2));
            idx = 1;
            for trial_n = 1:size(data_matrix, 1)
                for angle_n = 1:size(data_matrix, 2)
                    temp = zeros(1, size(data_matrix, 2));
                    for i = 1:size(data_matrix, 2)
                        temp(1, i) = mean(abs(templates(i, :) - (sum(squeeze(data_matrix(trial_n, angle_n, neurons_selected, start:stop)), 2)./(stop-start+1))'));
                    end
                    av_diff(idx, :) = temp;
                    idx = idx + 1;
                end
            end
        end
        
        function specific_neurons = mostSpecific(~, trial, start, stop, neuronsXangle)
            n_a = size(trial, 2);
            n_n = size(trial, 3);
            
            avg_activity = zeros(n_n, n_a);
            
            avg_activity = squeeze(mean(mean( ...
                           trial(:, :, :, start:stop) ...
                           , 4), 1))';
            sum_activity = sum(avg_activity, 2);
            avg_activity = avg_activity * 2;
            
            specific = zeros(n_n, n_a);
            specific = 2*(avg_activity - repmat(sum_activity, 1, n_a))./avg_activity;
            
            [~, specific_neurons] = sort(specific, 'descend');
            specific_neurons = specific_neurons(1:neuronsXangle, :);
        end
        
        % distribution of activity across neurons
        function distribution = firingDistribution(~, trial, start, stop)
            n_n = size(trial, 3);
            
            avg_activity = squeeze(mean(mean( ...
                           trial(:, :, :, start:stop) ...
                           , 4), 1))';
            
            sum_activity = sum(avg_activity, 1);
            
            distribution = avg_activity ./ ...
                           repmat(sum_activity, n_n, 1);
            
        end
        
        % sorted distribution of activity across angles
        function [pref_mag, pref_neuron, sum_activity] = anglePreference(~, trial, start, stop)
            n_a = size(trial, 2);
            
            avg_activity = squeeze(mean(mean( ...
                           trial(:, :, :, start:stop) ...
                           , 4), 1))';
            
            sum_activity = sum(avg_activity, 2);
            
            distribution = avg_activity ./ ...
                           repmat(sum_activity, 1, n_a);
                       
            [pref_mag, pref_neuron] = sort(distribution, 'descend');
            
        end
        
        function activity = overallActivity(~, trial, start, stop)
            activity = squeeze(sum(mean( ...
                       trial(:,:,:,start:stop) ...
                       , 1), 3));
        end
        
        
    end
    
    
end