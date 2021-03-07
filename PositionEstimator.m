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
        
        function W = calculateW(labels, A)
            c1 = 0;
            c2 = 0;
            for k = 2:size(labels, 2)
                c1 = c1 + (labels(:,k)*labels(:,k)');
                c2 = c2 + (labels(:,k-1)*labels(:,k)');
            end
            W = (1/size(labels, 2)-1)*(c1 - A*c2);
        end
        
        function H = calculateH()
            
        end
        
        function Q = calculateQ(samples, labels, H)
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