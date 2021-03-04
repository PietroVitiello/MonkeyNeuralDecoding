classdef Processing
    
    properties
        
    end
    
    methods
        function active_neurons = mostActive(trial, n, lower_bound, upper_bound)
            nargin
            
            if nargin < 4
                lower_bound = 1;
            end
            if nargin < 3
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
                        row = 1;
                    else
                        row = row + 1;
                    end
                    
                end   
            end 
        end
        
        
        
    end
    
end