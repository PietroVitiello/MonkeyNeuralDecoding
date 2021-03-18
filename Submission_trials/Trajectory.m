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
        
        function obj = objective_positions(~, hand)
            n_a = size(hand, 2);
            obj = zeros(2, n_a);
            obj = squeeze(mean(mean(...
                  hand(:, :, :, end-50:end), 4), 1))';
        end
        
        function [avgT, stdT] = averageTrajectory(~, hand)
            n_a = size(hand, 2);
            time_points = 320 : 20 : 571;
            avgT = zeros(n_a, 2, length(time_points));
            stdT = zeros(n_a, 2, length(time_points));
            
            avgT = squeeze(mean(hand(:,:,:,time_points), 1));
            stdT = squeeze(std(hand(:,:,:,time_points), 1));
            
        end
        
    end
end