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
        
        function avgT = averageTrajectory(~, trial)
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
                        eeg(:, i) = sum(trial(tr, a).spikes(:,bin_starts(i)-lag:bin_starts(i)+bin_size-1-lag)...
                                    .* repmat(((1.25).^(-(1:bin_size))), 98, 1), 2);
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
        
        
        
        
    end
end