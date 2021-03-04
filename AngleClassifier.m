classdef AngleClassifier
    
    properties
        
    end
    
    
    
    methods
        
        function classified_data = k_nn(~, training_data, labels, test_data, k)
            classified_data = zeros(1, length(test_data));
            [nearest_neighbours, dist] = knnsearch(training_data, test_data, 'K', k);
            for i = 1:length(test_data)
                nn_indices = nearest_neighbours(i, :);
                nn_labels = labels([nn_indices]);
                classified_data(i) = mode(nn_labels);
            end 
        end
        
        
        function [estimated_angles, true_angles] = likelihood(~, train_matrix, test_matrix)
            
            % training
            [~, n_angles, n_neurons, ~] = size(train_matrix);
            
            train_data = sum(train_matrix, 4);
            par = zeros(n_neurons, n_angles, 2);
            
            for neuron = 1:n_neurons
                for angle = 1:n_angles
                    [par(neuron, angle, :), ~] = mle(train_data(:, angle, neuron));
                end
            end
            
            % testing
            n_test_trials = size(test_matrix, 1)*n_angles;
            
            test_data = sum(test_matrix, 4);
            test_data = reshape(test_data, [n_test_trials, n_neurons]);
            
            estimated_angles = zeros(n_test_trials, 1);
            true_angles = zeros(n_test_trials, 1);
            
            like_a = zeros(n_test_trials, n_neurons, n_angles);
            
            for i = 1:n_test_trials
                for n = 1:n_neurons
                    for a = 1:n_angles
                        like_a(i, n, a) = normpdf(test_data(i, n), par(n,a,1), apr(n,a,2));
                    end
                end
                likelihood = prod(like_a(i,:,:));
                [~, estimated_angles(i)] = max(likelihood);
                true_angles(i) = floor((i-1) / size(test_matrix, 1));
            end
        end
        
        
    end
    
end