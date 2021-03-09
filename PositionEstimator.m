classdef PositionEstimator
    
    properties
        
    end
    
    methods
        
        function [x_estim, P_estim] = update(~, A, x_prev, H, Q, R, P_prev, obs)
            x_pred = A*x_prev;
            P_pred = A*P_prev*A' + Q;
            
            K_gain = P_pred*H'*(inv(H*P_pred*H' + R));
            x_estim = x_pred + K_gain*(obs - H*x_pred);
            P_estim = (eye(size(x_prev, 1), size(x_prev, 1)) - K_gain*H)*P_pred;
        end
        
        function [x_train, x_test] = getLabels(~, trial, delta, percent, win_size, deriv, start)
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
            x_train: 2*(deriv+1) x (total time steps) labels signal for training
            x_test: 2*(deriv+1) x (total time steps) labels for testing
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
            if nargin < 7
                start = 1;
            end
            if nargin < 6
                deriv = 1;
            end
            if nargin < 5
                win_size = 1;
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
        
        function [state0, eeg_train, eeg_test, x_train, x_test] = getDataset(~, trial, lag, percent, start)
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
            eeg_train: stimulus signal for training
            eeg_test: stimulus signal for testing
            x_train: labels signal for training
            x_test: labels for testing
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
            if nargin < 5
                start = 301;
            end
            
            first_t = start - lag;
            [n_tr, n_a] = size(trial); % #trials, #angles
            
            rand_id = randperm(n_tr);
            n_train = floor(percent * n_tr / 100);
            n_test = n_tr-n_train;
            
            state0 = zeros(2, n_a);
            eeg_train = cell(n_a, n_train, 1);
            x_train = cell(n_a, n_train, 1);
            eeg_test = cell(n_a, n_test, 1);
            x_test = cell(n_a, n_test, 1);
            
            for a = 1:n_a
                for i_tr = 1:n_tr
                    tr = rand_id(i_tr);
                    
                    average_handDisp = mean(diff(trial(tr, a).handPos(1:2,1:first_t-1)));
                    state0(:, a) = state0(:, a) + average_handDisp/n_a;
                    
                    eeg = trial(tr, a).spikes(:,first_t:end-lag);
                    x = diff(trial(tr, a).handPos(1:2,start:end));
                    
                    if i_tr <= n_train
                        eeg_train{a,i_tr,1} = eeg;
                        x_train{a,i_tr,1} = x;
                    else
                        eeg_test{a,i_tr-n_train,1} = eeg;
                        x_test{a,i_tr-n_train,1} = x;
                    end

                    eeg = [];
                    x = [];
                    
                end
            end    
        end
        
        function A = calculateA(~, x_cell)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            The purpose of this function is to create the dynamics matrix
            for the labels
            
            -input
            x: (label dimensions) x (time steps)
            
            -output
            A: labels dynamics matrix
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
            d = size(x_cell{1,1}, 1);
            sum1 = zeros(d);
            sum2 = zeros(d);
            
            for tr = 1:size(x_cell, 2)
                x = x_cell{1,tr};
                M = size(x, 2);
                for k = 2:M
                    sum1 = sum1 + (x(:,k) * x(:,k-1)');
                    sum2 = sum2 + (x(:,k-1) * x(:,k-1)');
                end
            end
            
            A = sum1/sum2;
        end
        
        function W = calculateW(~, labels, A)
            c1 = zeros(size(labels, 1));
            c2 = zeros(size(labels, 1));
            for k = 2:size(labels, 2)
                c1 = c1 + (labels(:,k)*labels(:,k)');
                c2 = c2 + (labels(:,k-1)*labels(:,k)');
            end
            W = (1/(size(labels, 2)-1))*(c1 - A*c2);
        end
        
        function H = calculateH(~, z, x)
            M = size(x, 2);
            sum1 = zeros(size(z, 1), size(x, 1));
            sum2 = zeros(size(x,1));
            
            for k = 1:M
                sum1 = sum1 + (z(:, k)*x(:, k)');
                sum2 = sum2 + (x(:, k)*x(:, k)');
            end
            
            H = sum1/sum2;
        end
        
        function Q = calculateQ(~, samples, labels, H)
            c1 = zeros(size(samples, 1));
            c2 = zeros(size(labels, 1), size(samples, 1));
            for k = 2:size(labels, 2)
                c1 = c1 + (samples(:,k)*samples(:,k)');
                c2 = c2 + (labels(:,k)*samples(:,k)');
            end
            Q = (1/size(labels, 2))*(c1 - H*c2);
        end       
        
        function [A, W, H, Q] = computeDynamics(obj, x, z)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            The purpose of this function is to return the dynamics and
            covariance matrices of the system
            
            -input
            x: (number of angles)x1 cell with each cell being
               (label dimensions) x (time steps)
            z: (number of angles)x1 cell with each cell being
               (number of neurons) x (time steps)
            
            -output
            A: labels dynamics matrix
            W: labels noise covariance
            H: stimulus dynamics matrix
            A: stimulus noise covariance
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
            A = [];
            W = [];
            H = [];
            Q = [];
            for a = 1:size(x)
                A = cat(3, A, obj.calculateA(x(a,:)));
                W = cat(3, W, obj.calculateW(x{a}, A(:,:,end)));
                H = cat(3, H, obj.calculateH(z{a}, x{a}));
                Q = cat(3, Q, obj.calculateQ(z{a}, x{a}, H(:,:,end)));
            end
        end
        
    end
end