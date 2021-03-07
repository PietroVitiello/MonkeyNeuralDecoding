classdef PositionEstimator
    
    properties
        
    end
    
    
    
    methods
        
        function A = calculateA(x, M)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            The purpose of this function is to create a list of neurons
            that have a high activity in the first 320ms to use for the 
            angle classifier
            
            -input
            x: trials x neurons
            n: number of most active neurons taken for each angle without
               repetition
            lower_bound: lower bound of the sample window to consider
            upper_bound: upper bound of the sample window to consider
            
            -output
            active_neurons: 1 x n*(number of angles) list of neurons
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
            l = size(x, 1);
            sum = zeros(l);
            
            for k = 2:M
                sum = sum + (x)
            end
            
        end
        
        function W = calculateW()
            
        end
        
        function H = calculateH()
            
        end
        
        function Q = calculateQ()
            
        end
        
    end
    
end