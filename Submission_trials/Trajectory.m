classdef Trajectory
    
    properties
        
    end
    
    methods
        
        function lens = movementDuration(~, trial)
            n_a = size(trial, 2);
            n_tr = size(trial, 1);
            lens = zeros(1, n_a);
            for a = 1:n_a
                for tr = 1:n_tr
                    lens(a) = lens(a) + size(trial(tr,a).spikes - 400,2)/n_tr;
                end
            end
        end
        
        function objective_positions(~, hand)
            n_a = size(hand, 2);
            obj = zeros(2, n_a);
            for a = 1:n_a
                obj(:, a) = mean(mean(...
                            hand(:, a, :, end-50:end), 4), 1);
            end
        end
        
        function [avgT, stdT] = averageTrajectory(~, trial)
            n_a = size(hand, 2);
            time_points = 300 : 20 : 571;
            avgT = zeros(n_a, 2, length(time_points));
            stdT = zeros(n_a, 2, length(time_points));
            
            for a = 1:n_a
                avgT(a, :, :) = mean(trial(:,a,:,time_points), 1);
                stdT(a, :, :) = std(trial(:,a,:,time_points), 1);
            end
        end
        
        
        
        
    end
end