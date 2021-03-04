classdef Processing
    
    properties
        trial = []
    end
    
    methods
        function active_neurons = mostActive(~, trial, n, lower_bound, upper_bound)
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
            
            if nargin < 5
                lower_bound = 1;
            end
            if nargin < 4
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
            
            active_neurons = [];
            col = 1;
            row = 1;
            
            while length(active_neurons) < n*size(all_active_neurons, 2)
                temp = all_active_neurons(row, col);
                
                if isempty(find(active_neurons == temp, 1))
                    active_neurons = [active_neurons temp];
                    
                    if length(active_neurons) == col*n
                        col = col + 1;
                        row = 0;
                    end
                end
                
                row = row + 1;
                
            end 
        end
        
        
        
        function [samples, labels] = create_dataset(trial, active_neurons)
            dataset = zeros(size(trial,1)*size(trial,2), length(active_neurons)+1);
            length_premotor = 320;
            traj_count = 0;
            for angle_n = 1:size(trial,2)
                for trial_n = 1:size(trial, 1)
                    traj_count = traj_count + 1;
                    temp = zeros(1, length(active_neurons));
                    for neuron_n = 1:length(active_neurons)
                        temp(neuron_n) = sum(trial(trial_n, angle_n).spikes(active_neurons(neuron_n), 1:length_premotor));
                        temp(length(active_neurons)+1) = angle_n;
                    end
                dataset(traj_count, :) = temp;
                end
            end
            rand_d = dataset(randperm(size(dataset, 1)), :);
            samples = rand_d(:, 1:length(active_neurons));
            labels = rand_d(:, length(active_neurons)+1);
        end
        
        
        
        
    end
    
end