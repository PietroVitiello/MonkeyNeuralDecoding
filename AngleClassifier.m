classdef AngleClassifier
    
    properties
        
    end
    
    methods
   
        function classified_data = k_nn(~, training_data, labels, test_data, k)
            classified_data = zeros(1, length(test_data));
            [nearest_neighbours, ~] = knnsearch(training_data, test_data, 'K', k);
            for i = 1:length(test_data)
                nn_indices = nearest_neighbours(i, :);
                nn_labels = labels([nn_indices]);
                classified_data(i) = mode(nn_labels);
            end 
        end
        
        function [train_metrics, test_metrics] = knn_classifier(~, train_samples, train_labels, test_samples, test_labels, k)
            Mdl = fitcknn(train_samples,train_labels,'NumNeighbors',k);
            rloss = resubLoss(Mdl);
            CVMdl = crossval(Mdl);
            kloss = kfoldLoss(CVMdl);
            train_metrics = [rloss, kloss];
            
            predictions = zeros(1, size(test_samples,1)); 
            for i = 1:length(test_samples)
                pred = predict(Mdl, test_samples(i, :));
                if pred == test_labels(i, :)
                    predictions(i) = 1;
                else
                    predictions(i) = 0;
                end
            end 
            
            test_metrics = (sum(predictions))/length(predictions);
        end
        
        
        function [estimated_angles, true_angles] = likelihood(~, train_matrix, test_matrix)
            
            %%% FITTING %%%
            [~, n_angles, n_neurons, ~] = size(train_matrix);
            
            train_data = sum(train_matrix, 4);
            par = zeros(n_neurons, n_angles, 2);
            
            for neuron = 1:n_neurons
                for angle = 1:n_angles
                    par(neuron, angle, :) = mle(train_data(:, angle, neuron));
                end
            end
            
            %%% TESTING %%%
            n_test_trials = size(test_matrix, 1)*n_angles;
            
            test_data = sum(test_matrix, 4);
            test_data = reshape(test_data, [n_test_trials, n_neurons]);
            
            estimated_angles = zeros(n_test_trials, 1);
            true_angles = zeros(n_test_trials, 1);
            
            like_a = zeros(n_test_trials, n_neurons, n_angles);
            
            for i = 1:n_test_trials
                for n = 1:n_neurons
                    for a = 1:n_angles
                        like_a(i, n, a) = normpdf(test_data(i, n), par(n,a,1), par(n,a,2));
                    end
                end
                likelihood = prod(like_a(i,:,:));
                [~, estimated_angles(i)] = max(likelihood);
                true_angles(i) = floor((i-1) / size(test_matrix, 1)) + 1;
            end
        end
        
        
        function [estimated_angles, true_angles] = multidimensional_mle(~, train_matrix, test_matrix)
            
            processor = Processing();
            
            %%% FITTING %%%
            d = size(train_matrix, 3);
            n_angles = size(train_matrix, 2);
            
            train_data = sum(train_matrix, 4);
            train_data = processor.switchDim(train_data);
            
            mu_ = mean(train_data, 1);
            size(mu_);

            sigma_inv_ = zeros(d, d, n_angles);
            scalar_ = zeros(n_angles, 1);
            for angle = 1:n_angles
                sigma = cov(train_data(:,:,angle));
                %add error message if not invertible?
                sigma_inv_(:,:,angle) = inv(sigma);
                
                scalar_(angle) = 1/sqrt(((2*pi)^d)*det(sigma));
            end
            
            %%% TESTING %%%
            n_test_trials = size(test_matrix, 1)*n_angles;
            
            test_data = sum(test_matrix, 4);
            test_data = reshape(test_data, [n_test_trials, d]);
            
            likelihood = zeros(n_angles, 1);
            estimated_angles = zeros(n_test_trials, 1);
            true_angles = zeros(n_test_trials, 1);
            
            for i = 1:n_test_trials
                for a = 1:n_angles
                    
                    x = test_data(i, :);
%                     x = mu_(1,:,a);
                    
                    mu = mu_(1,:,a);
                    sigma_inv = sigma_inv_(:,:,a);
                    scalar = scalar_(a);
                    
                    e = exp((-1/2)*dot((x-mu)*sigma_inv, (x-mu)'));
                    
                    likelihood(a) = scalar * e;
                    
                end
                [~, estimated_angles(i)] = max(likelihood);
                true_angles(i) = floor((i-1) / size(test_matrix, 1)) + 1;
            end
        end
        
        
    end
    
end