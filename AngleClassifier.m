classdef AngleClassifier
    
    properties
        
    end
    
    
    
    methods
        function classified_data = k_nn(~, training_data, test_data, k)
            classified_data = zeros(1, length(test_data));
            [nearest_neighbours, dist] = knnsearch(training_data, test_data, 'K', k);
            for i = 1:length(test_data)
                for j = 1:k
                    
                end
            end 
        end
    end
    
end