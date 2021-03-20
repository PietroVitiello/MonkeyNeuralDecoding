classdef Processing
    
    properties
      
    end
    
    methods

        function trial = clean_dataset(~, trial, silent_neurons)
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        The purpose of this function is to create a data matrix of structs 
        containing only neurons with a substantial firing rate across
        trials and angles
        -input
        trial: data as a matrix of structs
        silent_neurons: list of 9 most silent neurons
        -output
        trial: data as a matrix of structs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
            for angle_n = 1:size(trial, 2)
                for trial_n = 1:size(trial, 1)
                    for neuron_n = 1:length(silent_neurons)
                        trial(trial_n, angle_n).spikes(silent_neurons(neuron_n), :) = [];
                        %Spike trains of silent neurons are removed from
                        %the trial struct
                    end
                end
            end
        end
        
        
        
        
        function active_neurons = mostActive(~, trial, n, lower_bound, upper_bound, mode)
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
            
            if nargin < 4
                lower_bound = 1;
                mode = 'vector';
            end
            if nargin < 5
                upper_bound = 320;
                mode = 'vector';
            end

            window = lower_bound : upper_bound;
            %window contains the time instants prior to the onset
            %of movement (hence premotor)

            average_spike_trains = zeros(size(trial(1,1).spikes, 1), size(trial, 2));
            %average_spike_trains is a matrix containing time-averaged
            %firing rates summed over trials, each column represents one
            %angle and each row represents a neuron

            for angle_n = 1:size(trial, 2)
                for i = 1:size(trial, 1)
                    average_spike_trains(:,angle_n) = average_spike_trains(:,angle_n) + mean(trial(i, angle_n).spikes(:, 1:window), 2);
                    %The mean is computed across time in the premotor
                    %window and is then summed across trials, for each
                    %neuron
                end
            end

            [~, all_active_neurons] = sort(average_spike_trains, 'descend');
            %all_active_neurons is a matrix, each column represents one
            %angle and the neurons are ordered from the highest to lowest
            %time-averaged firing rate
            
            if (nargin > 4) && (strcmp(mode, 'matrix'))
                active_neurons = all_active_neurons(1:n, :);
            
            else
                active_neurons = [];
                col = 1;
                row = 1;

                while length(active_neurons) < n*size(all_active_neurons, 2)
                %active_neurons must be n*8 elements long 
                    temp = all_active_neurons(row, col);

                    if isempty(find(active_neurons == temp, 1))
                        active_neurons = [active_neurons temp];
                        %If the temp neuron is not contained in the 
                        %active_neurons list, it is added to it

                        if length(active_neurons) == col*n
                            col = col + 1;
                            row = 0;
                            %The ordered firing rates matrix is searched row by
                            %row (i.e. the whole top row -> the whole second
                            %row etc.)
                        end
                    end
                    row = row + 1;
                end
            end
        end
        
        
        function [samples, labels] = create_dataset(~, trial, active_neurons)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            The purpose of this function is to create a dataset of neurons
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
            dataset = zeros(size(trial,1)*size(trial,2), length(active_neurons)+1);
            %dataset is a matrix, the rows represent trials and the columns
            %represent the active neurons
            length_premotor = 320;
            traj_count = 0;
            for angle_n = 1:size(trial,2)
                for trial_n = 1:size(trial, 1)
                    traj_count = traj_count + 1;
                    temp = zeros(1, length(active_neurons));
                    for neuron_n = 1:length(active_neurons)
                        temp(neuron_n) = sum(trial(trial_n, angle_n).spikes(active_neurons(neuron_n), 1:length_premotor));
                        %Spikes are summed for each angle, trial and neuron
                        %in the first 320ms of each spike train ("premotor"
                        %window
                        temp(length(active_neurons)+1) = angle_n;
                        %The angle corresponding to a given trial is also
                        %recorded and added at the end of the row
                    end
                dataset(traj_count, :) = temp;
                end
            end
            rand_d = dataset(randperm(size(dataset, 1)), :);
            %The rows of the dataset are randomised to alternate trials of
            %different angles
            samples = rand_d(:, 1:length(active_neurons));
            %samples is a matrix containing spikes summed across time, the
            %rows represent trials and the columns represent active neurons 
            labels = rand_d(:, length(active_neurons)+1);
            %labels is a matrix containing the angles of reference for each
            %trial, the rows represent the trials
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
        

        function [firing_rates, av_labels] =  get_average_data(~,eeg_train, x_train)
            max_lengths = zeros(size(eeg_train,1), 1);
            for angle_n = 1:size(eeg_train, 1)
                max_ = -Inf;
                for trial_n = 1:size(eeg_train, 2)
                    if size(eeg_train{angle_n, trial_n},2) > max_
                        max_ = size(eeg_train{angle_n, trial_n},2);
                    end
                end
                max_lengths(angle_n) = max_;
            end
     
            firing_rates = cell(size(eeg_train,1), 1);
            av_labels = cell(size(x_train,1), 1);
            for angle_n = 1:size(eeg_train,1)
                fr_angle = zeros(size(eeg_train{angle_n, 1}, 1), max_lengths(angle_n));
                labels_angle = zeros(size(x_train{angle_n, 1}, 1), max_lengths(angle_n));
                for timebin = 1:max_lengths(angle_n)
                    fr_temp = zeros(size(eeg_train{angle_n, 1}, 1), 1);
                    fr_dividend = size(eeg_train,2);
                    if timebin <= max_lengths(angle_n)-1
                        label_temp = zeros(size(x_train{angle_n, 1}, 1), 1);
                        label_dividend = size(x_train, 2);
                    end
                    for trial_n = 1:size(eeg_train, 2)
                        if size(eeg_train{angle_n, trial_n}, 2) >= timebin
                            fr_temp = fr_temp + eeg_train{angle_n, trial_n}(:, timebin);
                        else
                            fr_dividend = fr_dividend - 1;
                        end
                        if size(x_train{angle_n, trial_n}, 2) >= timebin
                            label_temp = label_temp + x_train{angle_n, trial_n}(:, timebin);
                        else
                            label_dividend = label_dividend -1;
                        end
                    end
                    fr_angle(:, timebin) = fr_temp./fr_dividend;
                    if timebin <= max_lengths(angle_n)-1
                        labels_angle(:, timebin) = label_temp./label_dividend;
                    end
                end
                firing_rates{angle_n} = fr_angle;
                av_labels{angle_n} = labels_angle;
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
        
    end
end