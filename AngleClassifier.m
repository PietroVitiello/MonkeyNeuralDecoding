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
        
    end
    
end