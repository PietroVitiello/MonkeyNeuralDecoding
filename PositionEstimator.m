classdef PositionEstimator
    
    properties
        
    end
    
    
    
    methods
        
        function [eeg_train, eeg_test, x_train, x_test] = getLabels(~, trial, delta, percent, start)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            The purpose of this function is to create the training and
            testing datasets both for stimulus and labels. The datasets are
            created without separating between trials and with random
            permutation of the millisecond recordings
            
            -input
            trial: the given struct
            delta: time lag between stimulus and label in ms
            percent: percentage of training data
            start: to which sample start (optional)
            
            -output
            eeg_train: stimulus signal for trianing
            eeg_test: stimulus signal for testing
            x_train: labels signal for trianing
            x_test: labels for testing
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}

            if nargin < 5
                start = 1;
            end
            
            first_t = start + delta;
            [n_tr, n_a] = size(trial); % #trials, #angles
            
            eeg_train = cell(n_a, 1);
            x_train = cell(n_a, 1);
            eeg_test = cell(n_a, 1);
            x_test = cell(n_a, 1);
            
            eeg = [];
            x = [];
            for a = 1:n_a
                for tr = 1:n_tr
                    for t = first_t:size(trial(tr,a).spikes, 2)
                        eeg = [eeg trial(tr, a).spikes(:,t-delta)];
                        hand_disp = trial(tr, a).handPos(1:2,t)-trial(tr, a).handPos(1:2,t-1);
                        x = [x hand_disp];
                    end
                end
                
                n_train = floor(percent * size(eeg,2) / 100);
                rand_id = randperm(size(eeg,2));
                
                eeg_train{a,1} = eeg(:, rand_id(1:n_train));
                eeg_test{a,1} = eeg(:, rand_id(n_train+1:end));
                x_train{a,1} = x(:, rand_id(1:n_train));
                x_test{a,1} = x(:, rand_id(n_train+1:end));

                eeg = [];
                x = [];
            end    
        end
        
        function A = calculateA(~, x, M)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            The purpose of this function is to create the dynamics matrix
            for the labels
            
            -input
            x: (time steps) x (label dimesnions)
            
            -output
            A: labels dynamics matrix
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
            d = size(x, 1);
            sum1 = zeros(d);
            sum2 = zeros(d);
            
            for k = 2:M
                sum1 = sum1 + (x(:,k) * x(:,k-1)');
                sum2 = sum2 + (x(:,k-1) * x(:,k-1)');
            end
            
            A = sum1/sum2;
        end
        
        function W = calculateW(~, labels, A)
            c1 = 0;
            c2 = 0;
            for k = 2:size(labels, 2)
                c1 = c1 + (labels(:,k)*labels(:,k)');
                c2 = c2 + (labels(:,k-1)*labels(:,k)');
            end
            W = (1/size(labels, 2)-1)*(c1 - A*c2);
        end
        
        function H = calculateH(~, z, x, M)
            sum1 = zeros(size(z, 1), size(x, 1));
            sum2 = zeros(size(z, 1), size(x, 1));
            
            for k = 1:M
                sum1 = sum1 + (z(k, :)*x(k, :)');
                sum2 = sum2 + (x(k, :)*x(k, :)');
            end
            
            H = sum1/sum2;
        end
        
        function Q = calculateQ(~, samples, labels, H)
            c1 = 0;
            c2 = 0;
            for k = 2:size(labels, 2)
                c1 = c1 + (samples(:,k)*samples(:,k)');
                c2 = c2 + (labels(:,k)*samples(:,k)');
            end
            Q = (1/size(labels, 2))*(c1 - H*c2);
        end
        
    end
    
end