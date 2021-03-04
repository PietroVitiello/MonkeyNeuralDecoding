classdef Processing
    
    properties
        
    end
    
    methods
        function active_neurons = mostActive(trial, lower_bound, upper_bound)
            
            if nargin < 3
                lower_bound = 1;
            end
            if nargin < 2
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
            [~, active_neurons] = sort(average_spike_trains, 'descend');
        end
        
        function dataset = create_dataset(trial, active_neurons)
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
        end
    end
    
end