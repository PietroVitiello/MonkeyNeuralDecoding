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
        
        function Mdl = knn_classifier_(~, k, varargin)
            train_samples = varargin{1};
            train_labels = varargin{2};
            Mdl = fitcknn(train_samples,train_labels,'NumNeighbors',k);
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
            
            mu_ = squeeze(mean(train_data, 1));
            size(mu_);

            sigma_inv_ = zeros(d, d, n_angles);
            scalar_ = zeros(n_angles, 1);
            for angle = 1:n_angles
                sigma = cov(train_data(:,:,angle));
                %add error message if not invertible?
                eig(sigma);
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
                    
                    mu = mu_(:,a)';
                    sigma_inv = sigma_inv_(:,:,a);
                    scalar = scalar_(a);
                    
                    e = exp((-1/2)*dot((x-mu)*sigma_inv, (x-mu)'));
                    
                    likelihood(a) = scalar * e;
                    
                end
                [~, estimated_angles(i)] = max(likelihood);
                true_angles(i) = floor((i-1) / size(test_matrix, 1)) + 1;
            end
        end
        
        
        
        function [estimated_angles, true_angles] = multidimensional_mlev2(~, train_matrix, test_matrix, bin_size)
            
            [tr_train, n_angles, d, ~] = size(train_matrix);
            
%             ind = reshape(repmat(...
%                   [1:d/bin_size]', 1, bin_size)' ...
%                   , d, 1);
%               
%             permute(repmat(ind, 1, n_angles, tr_train), [3 2 1]);
            
            %%% FITTING %%%
            
            train_data = sum(train_matrix, 4);
            train_data = movsum(train_data, bin_size, 3, 'Endpoints','discard');
            train_data = train_data(:,:,1:bin_size:end);
            
            train_data = permute(train_data, [1 3 2]);
            
            mu_ = squeeze(mean(train_data, 1));
            size(mu_);

            sigma_inv_ = zeros(d/bin_size, d/bin_size, n_angles);
            scalar_ = zeros(n_angles, 1);
            for angle = 1:n_angles
                sigma = cov(train_data(:,:,angle));
                %add error message if not invertible?
                eig(sigma);
                sigma_inv_(:,:,angle) = inv(sigma);
                
                scalar_(angle) = 1/sqrt(((2*pi)^d)*det(sigma));
            end
            
            %%% TESTING %%%
            tr_test = size(test_matrix, 1);
            n_test_trials = tr_test*n_angles;
            
            test_data = sum(test_matrix, 4);
            test_data = movsum(test_data, bin_size, 3, 'Endpoints','discard');
            test_data = test_data(:,:,1:bin_size:end);
            
            test_data = reshape(test_data, [n_test_trials, d/bin_size]);
            
            likelihood = zeros(n_angles, 1);
            estimated_angles = zeros(n_test_trials, 1);
            true_angles = zeros(n_test_trials, 1);
            
            for i = 1:n_test_trials
                for a = 1:n_angles
                    
                    x = test_data(i, :);
%                     x = mu_(1,:,a);
                    
                    mu = mu_(:,a)';
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
        
        
        function templates = firingTemplate2(~, trial, stop, start, last_scale)
            
            n_a = size(trial, 2);
            n_n = size(trial, 3);
            scaling = [last_scale : (1-last_scale)/(stop-start) : 1]';
            
            trial2 = zeros(size(trial, 1), n_a, n_n);
            templates = zeros(n_a, n_n, 1);
            
            for tr = 1:size(trial, 1)
                for a = 1:n_a
                    for n = 1:n_n
                        trial(tr,a,n,start:stop);
                        trial2(tr,a,n) = sum(squeeze(trial(tr,a,n,start:stop)) .* scaling);
                    end
                end
            end
            
            templates = squeeze(mean(trial2, 1));
            
%             templates = squeeze(mean(sum(...
%                         trial(:,:,:,start:stop) .*...
%                         repmat(scaling, 98, )...
%                         , 4), 1));
                    
        end
        
        
        function angle = findSimilarAngle(~, templates, spikes, stop, start)
            
            B = mean(spikes(:, start:stop), 2)';
            [~, angle] = min(sum(abs(templates - ...
                         repmat(B, 8, 1))...
                         , 2));
        end
        
        % Worse, don't use
        function angle = findSimilarAngle1_2(~, templates, spikes, stop, start)
            
            B = mean(spikes(:, start:stop), 2)';
            [~, angle] = min(prod((abs(templates - ...
                         repmat(B, 8, 1) ...
                         )+1).*10, 2));
        end
        
        function angle = findSimilarAngle2(~, templates, spikes, stop, start, last_scale)
            
            scaling = last_scale : (1-last_scale)/(stop-start) : 1;
            B = sum( ...
                spikes(:,start:stop) .* repmat(scaling, size(spikes, 1), 1) ...
                , 2)';
            
            [~, angle] = min(sum(abs(templates - ...
                         repmat(B, 8, 1))...
                         , 2));
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
        
        function angle = aply_anglePreference(~, pref_mag, pref_neuron, sum_activity, spikes, start, stop)
           
            spikesAdistribution = repmat(mean(spikes(:, start:stop), 2) ...
                                  ./ sum_activity, 1, 8);
                              
            colI = repmat(1:8, size(pref_mag, 1), 1);
            chosen_mag = spikesAdistribution( ...
                         sub2ind(size(spikesAdistribution) ...
                         , pref_neuron, colI));
                     
            assignin('base', 'a', spikesAdistribution);
            assignin('base', 'b', pref_neuron);
            assignin('base', 'c', chosen_mag);
            
            [~, angle] = min(sum(abs(pref_mag - ...
                         chosen_mag)...
                         , 1), [], 2);
            
        end
        
        function angle = bin_activity(~, pref_mag, pref_neuron, sum_activity, spikes, start, stop)
           
            spikesAdistribution = repmat(mean(spikes(:, start:stop), 2) ...
                                  ./ sum_activity, 1, 8);
                              
            colI = repmat(1:8, size(pref_mag, 1), 1);
            chosen_mag = spikesAdistribution( ...
                         sub2ind(size(spikesAdistribution) ...
                         , pref_neuron, colI));
                     
            [~, angle] = max(sum( ...
                         chosen_mag)...
                         , [], 2);
            
        end
        
        function angle = binSimilarityDistributions(~, templates, angle_distributions, pref_neuron, sum_activity, spikes, start, stop)
            
            B = mean(spikes(:, start:stop), 2)';
            similar_firing = sum(abs(templates - ...
                             repmat(B, 8, 1))...
                             , 2);
            
            distribution = repmat( ...
                           B' / sum(B') ...
                           , 1, 8);
                       
            similar_dist = sum(abs(angle_distributions - ...
                           distribution)...
                           , 1)';
                       
            spikesAdistribution = repmat(B ...
                                  ./ sum_activity, 1, 8);
                              
            colI = repmat(1:8, size(pref_neuron, 1), 1);
            chosen_mag = spikesAdistribution( ...
                         sub2ind(size(spikesAdistribution) ...
                         , pref_neuron, colI));
                     
            binned_activity = sum( ...
                              chosen_mag);
                          
            [~, angle] = min(similar_firing + similar_dist - 0.005*binned_activity');
        end
        
        
        function par = neuronDistribution_mle(~, trial, start, stop)
            n_n = size(trial, 3);
            n_a = size(trial, 2);
            
            avg_activity = mean( ...
                           trial(:, :, :, start:stop) ...
                           , 4);
            
            sum_activity = sum(avg_activity, 3);
            
            distribution = avg_activity ./ ...
                           repmat(sum_activity, 1, 1, n_n);
            
            par = zeros(n_n, n_a, 2);
%             par(:,:,1) = arrayfun(@(aI, nI) mle(distribution(:,aI,nI)), ...
%                            1:n_a, 1:n_n);
            for n = 1:n_n
                for a = 1:n_a
                    par(n, a, :) = mle(distribution(:, a, n));
                end
            end
            
        end
        
        
        function angle = apply_nDistribution_mle(~, spikes, par, start, stop)
            
            n_n = size(par, 1);
            n_a = size(par, 2);
            
            avg_activity = mean( ...
                           spikes(:, start:stop) ...
                           , 2);
            
            sum_activity = sum(avg_activity);
            
            distribution = avg_activity / sum_activity;
            
            like_a = zeros(n_n, n_a);
            for n = 1:n_n
                for a = 1:n_a
                    like_a(n, a) = normpdf(distribution(n), par(n,a,1), par(n,a,2));
                end
            end
            silent_n = [8, 24, 25, 38, 42, 49, 52, 54, 73, 74, 76];
            like_a(silent_n,:) = [];
            like_a(isnan(like_a)) = 0.000001;
            likelihood = prod(like_a);
            [~, angle] = max(likelihood);
            
        end
        
        
        function [vector, angle] = visualizeSimilarityVectors(~, templates, spikes, stop, start)
            
            B = mean(spikes(:, start:stop), 2)';
            vector = sum(abs(templates - ...
                      repmat(B, 8, 1))...
                      , 2);
            [~, angle] = min(vector);
        end
        
        
        function angle = SimilarityMLE(~, templates, par, spikes, stop, start)
            n_n = size(par, 1);
            n_a = size(par, 2);
            
            avg_activity = mean( ...
                           spikes(:, start:stop) ...
                           , 2);
            
            sum_activity = sum(avg_activity);
            
            distribution = avg_activity / sum_activity;
            
            like_a = zeros(n_n, n_a);
            for n = 1:n_n
                for a = 1:n_a
                    like_a(n, a) = normpdf(distribution(n), par(n,a,1), par(n,a,2)) / 100;
%                     if like_a(n, a) < 1
%                         like_a(n, a) = 1;
%                     end
                end
            end
            silent_n = [8, 24, 25, 38, 42, 49, 52, 54, 73, 74, 76];
            like_a(silent_n,:) = [];
            like_a(isnan(like_a)) = 1;
            likelihood = prod(like_a);
            likelihood(1);
            
            vector = sum(abs(templates - ...
                      repmat(avg_activity', 8, 1))...
                      , 2);
            vector = (vector / sum(vector));
%             vector = 1 - vector
            
            [~, angle] = max(likelihood - vector');
            
            
        end
        
        function angle = compareAverage(~, training_average, test_trial, neurons)
            MSE_vector = zeros(1, 8);
            for i = 1:size(MSE_vector, 2)
                MSE_vector(1, i) = immse(training_average(:, :, i), test_trial(neurons, :));
            end
            [~, angle] = min(MSE_vector);
            MSE_vector;
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
        
        
        function templates = firingTemplate_2n3Dv1_2(~, trial, stop, start, bin_size)
            [~, n_a, n_n, ~] = size(trial);
            n_bin = (stop-start+1)/bin_size;
            trial = trial(:,:,:,start:stop);
            
            trial = squeeze(mean(trial, 1));
            
            templates = movsum(trial, bin_size, 3, 'Endpoints','discard');
            templates = reshape( ...
                        templates(:,:,1:bin_size:end) ...
                        , n_a, n_n*n_bin);
            
            templates = cat(2, mean(trial, 3), templates);
            templates = cat(2, ...
                        templates(:,1:n_n) ./ ...
                        repmat(sum(templates(:,1:n_n),2) ...
                        , 1,n_n)...
                        , templates);
                    
        end
        
        
        
        
        
        
    end
    
end