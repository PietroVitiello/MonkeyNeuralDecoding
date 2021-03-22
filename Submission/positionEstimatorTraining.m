function [modelParameters] = positionEstimatorTraining(training_data)
  % Arguments:
 
  % - training_data:
  %     training_data(n,k)              (n = trial id,  k = reaching angle)
  %     training_data(n,k).trialId      unique number of the trial
  %     training_data(n,k).spikes(i,t)  (i = neuron id, t = time)
  %     training_data(n,k).handPos(d,t) (d = dimension [1-3], t = time)
  
  % Return Value:
  
  % - modelParameters:
  %     - A
  %     - W
  %     - H
  %     - Q
  %     - classifier
  %     - initial parameters
  
  
  processor = Processing();
  trj = Trajectory();
  a_classifier = AngleClassifier();

    trials = zeros(size(training_data,1), 8, 98, 571);
    pos = zeros(size(training_data,1), 8, 2, 571);
    final_pos = zeros(size(training_data,1), 8, 2, 100);
    for trial_n = 1:size(training_data,1)
        for angle_n = 1: size(training_data,2)
            trials(trial_n, angle_n, :, :) = training_data(trial_n, angle_n).spikes(:, 1:571);
            pos(trial_n, angle_n, :, :) = training_data(trial_n, angle_n).handPos(1:2, 1:571);
            final_pos(trial_n,angle_n,:,:) = training_data(trial_n, angle_n).handPos(1:2, end-99:end);
        end
    end
    start = squeeze(mean(pos(:,:,:,1),1))';
    
    [n_tr, n_a, n_n, ~] = size(trials);
    n_bin = 300/100;
    trial = trials(:,:,:,1:300);

    templates2 = movsum(trial, 100, 4, 'Endpoints','discard');
    templates2 = reshape( ...
                templates2(:,:,:,1:100:end) ...
                , n_tr, n_a, n_n*n_bin);

    templates2 = cat(3, mean(trial, 4), templates2);
    templates2 = cat(3, ...
                templates2(:,:,1:n_n) ./ ...
                repmat(sum(templates2(:,:,1:n_n),3) ...
                , 1,1,n_n)...
                , templates2);

    templates2 = squeeze(mean(templates2, 1));
    
    n_bin = 400/100;
    trial = trials(:,:,:,1:400);

    templates4 = movsum(trial, 100, 4, 'Endpoints','discard');
    templates4 = reshape( ...
                templates4(:,:,:,1:100:end) ...
                , n_tr, n_a, n_n*n_bin);

    templates4 = cat(3, mean(trial, 4), templates4);
    templates4 = cat(3, ...
                templates4(:,:,1:n_n) ./ ...
                repmat(sum(templates4(:,:,1:n_n),3) ...
                , 1,1,n_n)...
                , templates4);

    templates4 = squeeze(mean(templates4, 1));
    
    n_bin = 500/125;
    trial = trials(:,:,:,1:500);

    templates6 = movsum(trial, 125, 4, 'Endpoints','discard');
    templates6 = reshape( ...
                templates6(:,:,:,1:125:end) ...
                , n_tr, n_a, n_n*n_bin);

    templates6 = cat(3, mean(trial, 4), templates6);
    templates6 = cat(3, ...
                templates6(:,:,1:n_n) ./ ...
                repmat(sum(templates6(:,:,1:n_n),3) ...
                , 1,1,n_n)...
                , templates6);

    templates6 = squeeze(mean(templates6, 1));
    
    time_points = 320 : 20 : 571;
    avgT = zeros(n_a, 2, length(time_points));
    stdT = zeros(n_a, 2, length(time_points));

    avgT = squeeze(mean(pos(:,:,:,time_points), 1));
    stdT = squeeze(std(pos(:,:,:,time_points), 1));
    
    obj = zeros(2, n_a);
    obj = squeeze(mean(mean(...
          final_pos, 4), 1))';

    modelParameters.initial = start;
    modelParameters.traces = avgT;
    modelParameters.deviation = stdT;
    modelParameters.objectives = obj;
    modelParameters.templates2 = templates2;
    modelParameters.templates4 = templates4;
    modelParameters.templates6 = templates6;
  
end