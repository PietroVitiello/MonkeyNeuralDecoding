classdef PositionEstimator
    
    properties
        
    end
    
    methods
        
        function [x_estim, P_estim] = update(~, A, x_prev, H, Q, R, P_prev, obs)
            x_pred = A*x_prev;
            P_pred = A*P_prev*A' + R;
            
            K_gain = P_pred*H'*(inv(H*P_pred*H' + Q));
            x_estim = x_pred + K_gain*(obs - H*x_pred);
            P_estim = (eye(size(x_prev, 1)) - K_gain*H)*P_pred;
        end
        
        
        
        
        
        
        
        
        function usable_data = apply_dataset_coala(~, data, lag, bin_size)
            
            [n_n, len] = size(data);
            first_t = len - 19 - lag;
            n_bin = floor(len / bin_size);
            bin_starts = first_t: bin_size: (first_t+(n_bin-1)*bin_size);
            
            usable_data = zeros(n_n, n_bin);
            for i = 1:n_bin
                usable_data(:, i) = mean(data(:,bin_starts(i):bin_starts(i)+bin_size-1),2);
            end
        end
        
        function [state0, eeg_train, eeg_test, x_train, x_test] = ferromagnetico(~, trial, lag, bin_size, order, percent, start)
            if nargin < 7
                start = 301;
            end
            
            first_t = start - lag;
            [n_tr, n_a] = size(trial); % #trials, #angles
            n_n = size(trial(1,1).spikes, 1); % #neurons
            
            rand_id = randperm(n_tr);
            n_train = floor(percent * n_tr / 100);
            n_test = n_tr-n_train;
            
            state0 = zeros(2*order, n_a);
            eeg_train = cell(n_a, n_train, 1);
            x_train = cell(n_a, n_train, 1);
            eeg_test = cell(n_a, n_test, 1);
            x_test = cell(n_a, n_test, 1);
            
            for a = 1:n_a
                for i_tr = 1:n_tr
                    tr = rand_id(i_tr);
                    len = size(trial(tr, a).spikes, 2);
                    n_bin = floor((len-start) / bin_size);
                    bin_starts = start: bin_size: (start+(n_bin-1)*bin_size);
                    
                    eeg = zeros(n_n, n_bin);
                    for i = 1:n_bin
                        eeg(:, i) = mean(trial(tr, a).spikes(:,bin_starts(i)-lag:bin_starts(i)+bin_size-1-lag),2);
                    end
                    
                    x = trial(tr, a).handPos(1:2,bin_starts+(bin_size - 1));
                    for o = 1 : order
                        x_temp = [];
                        for i = 1:n_bin
                            starting = bin_starts(i) - o;
                            x_temp = [x_temp, ...
                                      mean(diff(trial(tr, a).handPos(1:2,starting:...
                                      (bin_starts(i)+bin_size-1)), o, 2),2)];
                        end
                        
                        x = [x; x_temp];
                        
                        average_state = mean(diff(trial(tr, a).handPos(1:2,1:first_t-1),o,2),2);
                        state0((1:2)*o, a) = state0((1:2)*o, a) + average_state/n_tr;
                    end
                    
                    if i_tr <= n_train
                        eeg_train{a,i_tr,1} = eeg;
                        x_train{a,i_tr,1} = x;
                    else
                        eeg_test{a,i_tr-n_train,1} = eeg;
                        x_test{a,i_tr-n_train,1} = x;
                    end
                    
                end
            end
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
        
        function usable_data = apply_dataset(~, data, lag, bin_size, start)
            if nargin < 5
                start = 301;
            end
            
            first_t = start - lag;
            [n_n, len] = size(data);
            
            rows = n_n*2;
            usable_data = zeros(rows, len-start+1);
            usable_data(1:2:rows, :) = data(:,first_t:end-lag);
            for history = 1:bin_size
                usable_data(2:2:rows, :) = usable_data(2:2:rows, :) + ...
                                           data(:,(first_t-history):(end-lag-history)) / bin_size;
            end
        end
        
        function [state0, eeg_train, eeg_test, x_train, x_test] = sayonara(~, trial, lag, bin_size, order, percent, start)
            if nargin < 7
                start = 301;
            end
            
            first_t = start - lag;
            [n_tr, n_a] = size(trial); % #trials, #angles
            n_n = size(trial(1,1).spikes, 1); % #neurons
            
            rand_id = randperm(n_tr);
            n_train = floor(percent * n_tr / 100);
            n_test = n_tr-n_train;
            
            state0 = zeros(2*(order+1), n_a);
            eeg_train = cell(n_a, n_train, 1);
            x_train = cell(n_a, n_train, 1);
            eeg_test = cell(n_a, n_test, 1);
            x_test = cell(n_a, n_test, 1);
            
            for a = 1:n_a
                for i_tr = 1:n_tr
                    tr = rand_id(i_tr);
                    
                    len = size(trial(tr, a).spikes, 2);
                    rows = n_n*2;
                    eeg = zeros(rows, len-start+1);
                    eeg(1:2:rows, :) = trial(tr, a).spikes(:,first_t:end-lag);
                    for history = 1:bin_size
                        eeg(2:2:rows, :) = eeg(2:2:rows, :) + ...
                                           trial(tr, a).spikes(:,(first_t-history):(end-lag-history)) / bin_size;
                    end
                    
                    x = trial(tr, a).handPos(1:2,start:end);
                    for o = 1 : order
                        starting = start - o;
                        x = [x; ...
                            diff(trial(tr, a).handPos(1:2,starting:end), o, 2)];
                        
                        average_state = mean(diff(trial(tr, a).handPos(1:2,1:first_t-1),o+1,2),2);
                        state0((1:2)*(o+1), a) = state0((1:2)*(o+1), a) + average_state/n_tr;
                    end
                    
                    if i_tr <= n_train
                        eeg_train{a,i_tr,1} = eeg;
                        x_train{a,i_tr,1} = x;
                    else
                        eeg_test{a,i_tr-n_train,1} = eeg;
                        x_test{a,i_tr-n_train,1} = x;
                    end
                    
                end
            end
            state0(1:2, :) = [];
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
                    
                    average_handDisp = mean(diff(trial(tr, a).handPos(1:2,1:first_t-1),1,2),2);
                    state0(:, a) = state0(:, a) + average_handDisp/n_tr;
                    
                    eeg = trial(tr, a).spikes(1:20,first_t:end-lag);
                    x = diff(trial(tr, a).handPos(1:2,start:end), 1, 2);
                    
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
                
                sum1 = sum1 + x(:,2:M) * x(:,1:M-1)';
                sum2 = sum2 + x(:,1:M-1) * x(:,1:M-1)';
                
            end
            
            A = sum1/sum2;
        end
        
        function W = calculateW(~, x_cell, A)
            d = size(x_cell{1,1}, 1);
            
            c1 = zeros(d);
            c2 = zeros(d);
            MM = 0;
            for tr = 1:size(x_cell, 2)
                x = x_cell{1,tr};
                M = size(x, 2);
                MM = MM + M - 1;
                
                c1 = c1 + x(:,2:M)*x(:,2:M)';
                c2 = c2 + x(:,1:M-1)*x(:,2:M)';
                
                %W = W + (1/(M-1))*(c1 - A*c2)./size(x_cell, 2);
            end
            W = (1/MM)*(c1 - A*c2);
        end
        
        function H = calculateH(~, z_cell, x_cell)
            d_x = size(x_cell{1,1}, 1);
            d_z = size(z_cell{1,1}, 1);
            
            sum1 = zeros(d_z, d_x);
            sum2 = zeros(d_x);
            for tr = 1:size(x_cell, 2)
                x = x_cell{1,tr};
                z = z_cell{1,tr};
                M = size(x, 2); 
                
                sum1 = sum1 + z(:, 1:M)*x(:, 1:M)';
                sum2 = sum2 + x(:, 1:M)*x(:, 1:M)';
                
            end
            
            H = sum1/sum2;
        end
        
        function Q = calculateQ(~, z_cell, x_cell, H)
            d_x = size(x_cell{1,1}, 1);
            d_z = size(z_cell{1,1}, 1);
            
            c1 = zeros(d_z);
            c2 = zeros(d_x, d_z);
            Q = zeros(d_z);
            MM = 0;
            for tr = 1:size(x_cell, 2)
                x = x_cell{1,tr};
                z = z_cell{1,tr};
                M = size(x, 2);
                MM = MM + M;
                
                c1 = c1 + z(:,1:M)*z(:,1:M)';
                c2 = c2 + x(:,1:M)*z(:,1:M)';
                
%                 Q = Q + (1/M)*(c1 - H*c2)./size(x_cell, 2);
            end
            Q = (1/MM)*(c1 - H*c2);
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
                W = cat(3, W, obj.calculateW(x(a,:), A(:,:,end)));
                H = cat(3, H, obj.calculateH(z(a,:), x(a,:)));
                Q = cat(3, Q, obj.calculateQ(z(a,:), x(a,:), H(:,:,end)));
            end
        end
        
    end
end