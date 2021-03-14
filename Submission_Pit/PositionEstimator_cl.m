classdef PositionEstimator_cl
    
    properties
        
    end
    
    methods
        
        function [x_estim, P_estim] = update(~, A, x_prev, H, Q, R, P_prev, obs)
            x_pred = A*x_prev;
            P_pred = A*P_prev*A' + R;
            
            K_gain = P_pred*H'*(pinv(H*P_pred*H' + Q));
            x_estim = x_pred + K_gain*(obs - H*x_pred);
            P_estim = (eye(size(x_prev, 1)) - K_gain*H)*P_pred;
        end
        
        function usable_data = apply_ferromagnetico(~, data, lag, bin_size)
            
            [n_n, len] = size(data);
            first_t = len - 19 - lag;
            n_bin = floor((len-first_t) / bin_size);
            bin_starts = first_t: bin_size: (first_t+(n_bin-1)*bin_size);
            
            usable_data = zeros(n_n+1, n_bin);
            for i = 1:n_bin
                usable_data(2:end, i) = mean(data(:,bin_starts(i):bin_starts(i)+bin_size-1),2);
            end
            usable_data = [bin_starts; usable_data];
            
        end
        
        function [state0, eeg_train, eeg_test, x_train, x_test] = ferromagnetico(~, trial, lag, bin_size, order, percent, start)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Return state0 and the different datasets where the rows from the
            second to the last of eeg represent each one neuron and each
            column is the average firing rate in a ceratin bin (which
            depends on the bin_size). The labels are the position and any
            order of derivation of it, and they are as well the mean of the
            labels in the bins (apart from the position which has been
            taken as the last position of the bin).
            
            -input
            trial: the given struct
            delta: time lag between stimulus and label in ms
            bin_size: size of the bins (should be a factor of 20 otherwise
                      things might get messy)
            order: derivative order of the position
            percent: percentage of training data
            start: to which sample start (optional)
            
            -output
            state0: initial set of labels (apart from position)
            eeg_train: stimulus signal for training
            eeg_test: stimulus signal for testing
            x_train: labels signal for training
            x_test: labels for testing
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
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
                    
                    eeg = zeros(n_n+1, n_bin);
                    for i = 1:n_bin
                        eeg(2:end, i) = mean(trial(tr, a).spikes(:,bin_starts(i)-lag:bin_starts(i)+bin_size-1-lag),2);
                    end
                    eeg = [bin_starts; eeg];
                    
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
        
        function usable_data = apply_sayonara_2(~, data, lag, bin_size)
            
            [n_n, len] = size(data);
            first_t = len - 19 - lag;
            
            rows = n_n*2;
            usable_data = zeros(rows, 20);
%             usable_data = zeros(n_n, 20);
            usable_data(1:2:rows, :) = data(:,first_t:end-lag);
            for history = 1:bin_size
%                 usable_data(2:2:rows, :) = usable_data(2:2:rows, :) + ...
%                 usable_data = usable_data + ...
                usable_data(2:2:rows, :) = usable_data(2:2:rows, :) + ...
                                           (data(:,(first_t-history):(end-lag-history)) ...
                                           * 10 / 2^(history));
            end
            step = 1;
            time_steps = (first_t:1:first_t+19).*step;
            usable_data = [exp(time_steps) ; usable_data];
        end
        
        function [state0, eeg_train, eeg_test, x_train, x_test] = sayonara_2(~, trial, lag, bin_size, order, percent, start)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Return state0 and the different datasets. The odd rows of the
            eeg represent the neurons firing each millisecond, the even
            rows represent the average spiking of the last n milliseconds
            (where n is the bin_size) for each millisecond. The labels are 
            the position and any order of derivation of it, and they are
            as well the mean of the labels in the bins (apart from the
            position which has been taken as the last position of the bin).
            
            -input
            trial: the given struct
            delta: time lag between stimulus and label in ms
            bin_size: size of the bins (should be a factor of 20 otherwise
                      things might get messy)
            order: derivative order of the position
            percent: percentage of training data
            start: to which sample start (optional)
            
            -output
            state0: initial set of labels (apart from position)
            eeg_train: stimulus signal for training
            eeg_test: stimulus signal for testing
            x_train: labels signal for training
            x_test: labels for testing
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
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
%                     eeg = zeros(n_n, len-start+1);
                    eeg(1:2:rows, :) = trial(tr, a).spikes(:,first_t:end-lag);
                    for history = 1:bin_size
%                         eeg(2:2:rows, :) = eeg(2:2:rows, :) + ...
%                         eeg = eeg + ...
                        eeg(2:2:rows, :) = eeg(2:2:rows, :) + ...
                                           (trial(tr, a).spikes(:,(first_t-history):(end-lag-history)) ...
                                           * 10 / 2^(history)) ;
                    end
                    step = 1;
                    time_steps = (1:1:len-start+1) .* step;
                    eeg = [exp(time_steps) ; eeg];
                    
                    x = [];
                    for o = 0 : order
                        starting = start - (o + 1);
                        x = [x; ...
                            diff(trial(tr, a).handPos(1:2,starting:end), o+1, 2)];
                        
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

        function usable_data = apply_sayonara(~, data, lag, bin_size)
            
            [n_n, len] = size(data);
            first_t = len - 19 - lag;
            
            rows = n_n*2;
            usable_data = zeros(rows, 20);
%             usable_data = zeros(n_n, 20);
            usable_data(1:2:rows, :) = data(:,first_t:end-lag);
            for history = 1:bin_size
%                 usable_data(2:2:rows, :) = usable_data(2:2:rows, :) + ...
%                 usable_data = usable_data + ...
                usable_data(2:2:rows, :) = usable_data(2:2:rows, :) + ...
                                           (data(:,(first_t-history):(end-lag-history)) ...
                                           * 10 / bin_size);
            end
%             step = 1;
%             time_steps = (first_t:1:first_t+19).*step;
%             usable_data = [exp(time_steps) ; usable_data];
        end
        
        function [state0, eeg_train, eeg_test, x_train, x_test] = sayonara(~, trial, lag, bin_size, order, percent, start)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Return state0 and the different datasets. The odd rows of the
            eeg represent the neurons firing each millisecond, the even
            rows represent the average spiking of the last n milliseconds
            (where n is the bin_size) for each millisecond. The labels are 
            the position and any order of derivation of it, and they are
            as well the mean of the labels in the bins (apart from the
            position which has been taken as the last position of the bin).
            
            -input
            trial: the given struct
            delta: time lag between stimulus and label in ms
            bin_size: size of the bins (should be a factor of 20 otherwise
                      things might get messy)
            order: derivative order of the position
            percent: percentage of training data
            start: to which sample start (optional)
            
            -output
            state0: initial set of labels (apart from position)
            eeg_train: stimulus signal for training
            eeg_test: stimulus signal for testing
            x_train: labels signal for training
            x_test: labels for testing
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
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
%                     eeg = zeros(n_n, len-start+1);
                    eeg(1:2:rows, :) = trial(tr, a).spikes(:,first_t:end-lag);
                    for history = 1:bin_size
%                         eeg(2:2:rows, :) = eeg(2:2:rows, :) + ...
%                         eeg = eeg + ...
                        eeg(2:2:rows, :) = eeg(2:2:rows, :) + ...
                                           (trial(tr, a).spikes(:,(first_t-history):(end-lag-history)) ...
                                           * 10 / bin_size) ;
                    end
%                     step = 1;
%                     time_steps = (1:1:len-start+1) .* step;
%                     eeg = [exp(time_steps) ; eeg];
                    
                    x = [];
                    for o = 0 : order
                        starting = start - (o + 1);
                        x = [x; ...
                            diff(trial(tr, a).handPos(1:2,starting:end), o+1, 2)];
                        
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
            W = (1/MM).*(c1 - A*c2);
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
            Q = (1/MM).*(c1 - H*c2);
        end       
        
        function [A, W, H, Q] = computeDynamics_(obj, x, z)
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
        
        
        function [A, W, H, Q] = computeDynamics(~, x_cell, z_cell)
            
            A = [];
            W = [];
            H = [];
            Q = [];
            
            d_x = size(x_cell{1,1}, 1);
            d_z = size(z_cell{1,1}, 1);
            
            
            for a = 1:size(x_cell)
                
                A_sum1 = zeros(d_x);
                A_sum2 = zeros(d_x);
                W_sum1 = zeros(d_x);
                W_sum2 = zeros(d_x);
                H_sum1 = zeros(d_z, d_x);
                H_sum2 = zeros(d_x);
                Q_sum1 = zeros(d_z);
                Q_sum2 = zeros(d_x, d_z);
                
                W_temp = zeros(d_x);
                Q_temp = zeros(d_z);
                
                for tr = 1:size(x_cell, 2)
                    M = size(x_cell{a,tr}, 2);
                    
                    x1 = x_cell{a,tr}(:,1:M-1);
                    x2 = x_cell{a,tr}(:,2:M);
                    x11 = x_cell{a,tr}(:,1:M);
                    z11 = z_cell{a,tr}(:,1:M);

                    A_sum1 = A_sum1 + x2 * x1';
                    A_sum2 = A_sum2 + x1 * x1';
                    
                    H_sum1 = H_sum1 + z11 * x11';
                    H_sum2 = H_sum2 + x11 * x11';
                    
                end
                
                A = cat(3, A, A_sum1/A_sum2);
                H = cat(3, H, H_sum1/H_sum2);
                
                for tr = 1:size(x_cell, 2)
                    M = size(x_cell{a,tr}, 2);
                    
                    x1 = x_cell{a,tr}(:,1:M-1);
                    x2 = x_cell{a,tr}(:,2:M);
                    x11 = x_cell{a,tr}(:,1:M);
                    z11 = z_cell{a,tr}(:,1:M);
                    
                    W_sum1 = W_sum1 + x2 * x2';
                    W_sum2 = W_sum2 + x1 * x2';
                    
                    Q_sum1 = Q_sum1 + z11 * z11';
                    Q_sum2 = Q_sum2 + x11 * z11';
                    
                    W_temp = W_temp + ((1/M)*(W_sum1 - A(:,:,end)*W_sum2)./size(x_cell, 2));
                    Q_temp = Q_temp + ((1/M)*(Q_sum1 - H(:,:,end)*Q_sum2)./size(x_cell, 2));
                end
                W = cat(3, W, W_temp);
                Q = cat(3, Q, Q_temp);
                
            end
        end
        
        
        
        
    end
end