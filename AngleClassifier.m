classdef AngleClassifier
    
    properties
        
    end
    
    
    
    methods
        function classified_data = k_nn(~, training_data, labels, test_data, k)
            classified_data = zeros(1, length(test_data));
            [nearest_neighbours, dist] = knnsearch(training_data, test_data, 'K', k);
            for i = 1:length(test_data)
                nn_indices = nearest_neighbours(i);
                nn_labels = labels([nn_indices]);
                classified_data(i) = mode(labels_nn);
            end 
        end
    end
    
end