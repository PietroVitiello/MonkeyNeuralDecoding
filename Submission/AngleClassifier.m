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
        
        function [Mdl, train_metrics, test_metrics] = knn_classifier(~, k, varargin)
            % Arguments:
            %   - train_samples
            %   - train_labels
            %   - test_samples
            %   - test_labels
            
            train_samples = varargin{1};
            train_labels = varargin{2};
            Mdl = fitcknn(train_samples,train_labels,'NumNeighbors',k);
            rloss = resubLoss(Mdl);
            CVMdl = crossval(Mdl);
            kloss = kfoldLoss(CVMdl);
            train_metrics = [rloss, kloss];
            test_metrics = [];
            
            if length(varargin) > 3
                test_samples = varargin{3};
                test_labels = varargin{4};
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
        
        function templates = firingTemplate(~, trial, stop, start)
            
            n_a = size(trial, 2);
            n_n = size(trial, 3);
            templates = zeros(n_a, n_n, 1);
            
            templates = squeeze(mean(mean(...
                        trial(:,:,:,start:stop)...
                        , 4), 1));
                    
                    
        end
        
        function angle = findSimilarAngle(~, templates, spikes, stop, start)
            
            B = mean(spikes(:, start:stop), 2)';
            [~, angle] = min(sum(abs(templates - ...
                         repmat(B, 8, 1))...
                         , 2));
                     
        end
        
        % Worse, don't use
        function angle = findSimilarAngle2(~, templates, spikes, stop, start)
            
            B = mean(spikes(:, start:stop), 2)';
            [~, angle] = min(prod((abs(templates - ...
                         repmat(B, 8, 1) ...
                         )+1).*10, 2));
                     
        end
        
        function angle = likelyDistribution(~, angle_distributions, spikes, start, stop)
            
            avg_activity = mean( ...
                           spikes(:, start:stop) ...
                           , 2);
            
            distribution = repmat( ...
                           avg_activity / sum(avg_activity) ...
                           , 1, 8);
                       
            [~, angle] = min(sum(abs(angle_distributions - ...
                         distribution)...
                         , 1)');
            
        end
        
        function angle = similarityDistributions(~, templates, angle_distributions, spikes, start, stop)
            
            B = mean(spikes(:, start:stop), 2)';
            similar_firing = sum(abs(templates - ...
                             repmat(B, 8, 1))...
                             , 2);
                     
            avg_activity = mean( ...
                           spikes(:, start:stop) ...
                           , 2);
            
            distribution = repmat( ...
                           avg_activity / sum(avg_activity) ...
                           , 1, 8);
                       
            similar_dist = sum(abs(angle_distributions - ...
                           distribution)...
                           , 1)';
                     
            [~, angle] = min(similar_firing + similar_dist);
            
        end
        
           
        function angle = closestFiring(~, active_neurons, spikes, stop, start)
           
            spikes_t_average = mean(spikes(:, stop:start), 2);
            n = size(active_neurons, 1);
            [~, most_firing] = sort(spikes_t_average, 'descend');
            current_active = most_firing(1:n);
            
            similarity_matrix = zeros(size(active_neurons, 1), size(active_neurons, 2));
            for angle_n = 1:size(active_neurons, 2)
                similarity_matrix(:, angle_n) = active_neurons(:, angle_n) == current_active;
            end
            similarity_vector = sum(similarity_matrix, 1);
            [~, angle] = max(similarity_vector);
            
        end
        
        function templates = firingTemplate_3D(~, trial, stop, start, bin_size)
            trial = trial(:,:,:,start:stop);
            templates = movsum(trial, bin_size, 4, 'Endpoints','discard');
            templates = templates(:,:,:,1:bin_size:end);
            
            templates = squeeze(mean(templates, 1));
                    
        end
        
        function angle = findSimilarAngle_3D(~, templates, spikes, stop, start, bin_size)
            spikes = spikes(:,start:stop);
            B = movsum(spikes, bin_size, 2, 'Endpoints','discard');
            B = B(:,1:bin_size:end);
            
            angle_dissim = sum(abs(templates - ...
                           permute(repmat(B, 1, 1, 8),[3 1 2]) ...
                           ),2);
            bin_sum = sum(angle_dissim, 1);
            angle_dissim = angle_dissim ./ repmat(bin_sum,8,1);
            [~, angle] = min(sum(angle_dissim, 3));
            
        end
        
        function templates = firingTemplate_2n3D(~, trial, stop, start, bin_size)
            [n_tr, n_a, n_n, ~] = size(trial);
            n_bin = (stop-start+1)/bin_size;
            trial = trial(:,:,:,start:stop);
            
            templates = movsum(trial, bin_size, 4, 'Endpoints','discard');
            templates = reshape( ...
                        templates(:,:,:,1:bin_size:end) ...
                        , n_tr, n_a, n_n*n_bin);
            
            templates = cat(3, mean(trial, 4), templates);
            templates = cat(3, ...
                        templates(:,:,1:n_n) ./ ...
                        repmat(sum(templates(:,:,1:n_n),3) ...
                        , 1,1,n_n)...
                        , templates);
%             templates(:,:,n_n+1:end) = templates(:,:,n_n+1:end) ./ ...
%                                        repmat(templates(:,:,1:n_n), 1,1,n_bin);
            
            assignin('base', 'ww', repmat(templates(:,:,1:n_n), 1,1,n_bin))
                                   
            templates = squeeze(mean(templates, 1));
%             templates(isnan(templates)) = 0;
            size(templates);
            
            assignin('base', 'qq', templates')
                    
        end
        
        function angle = findSimilarAngle_2n3D(~, templates, spikes, stop, start, bin_size)
            n_n = size(spikes, 1);
            n_bin = (stop-start+1)/bin_size;
            spikes = spikes(:,start:stop);
            
            B = movsum(spikes, bin_size, 2, 'Endpoints','discard');
            B = reshape(B(:,1:bin_size:end), n_n*n_bin, 1);
            B = [mean(spikes, 2) ; B];
            B = [B(1:n_n) ./ repmat(sum( ...
                B(1:n_n)), n_n, 1); B];
%             size(B)
%             size(repmat(B(1:n_n), n_bin, 1))
%             B(n_n+1:end) = B(n_n+1:end) ./ repmat(B(1:n_n), n_bin, 1);
%             B(isnan(B)) = 0;
            
            [~, angle] = min(sum(abs(templates - ...
                         repmat(B, 1, 8)' ...
                         ), 2));
            
        end
        
    end
    
end