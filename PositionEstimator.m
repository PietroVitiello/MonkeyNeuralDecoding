classdef PositionEstimator
    
    properties
        
    end
    
    
    
    methods
        
        function A = calculateA(x, M)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            The purpose of this function is to create the dynamics matrix
            for the labels
            
            -input
            M: number of time points
            x: M x label_dimesnions
            
            -output
            A: dynamics matrix
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
            d = size(x, 1);
            sum1 = zeros(d);
            sum2 = zeros(d);
            
            for k = 2:M
                sum1 = sum1 + (x(k,:)' * x(k-1,:));
                sum2 = sum2 + (x(k-1,:)' * x(k-1,:));
            end
            
            A = sum1/sum2;
            
        end
        
        function W = calculateW()
            
        end
        
        function H = calculateH(z, x, M)
            sum1 = zeros(size(z, 1), size(x, 1));
            sum2 = zeros(size(z, 1), size(x, 1));
            
            for k = 1:M
                sum1 = sum1 + (z(k, :)*x(k, :)');
                sum2 = sum2 + (x(k, :)*x(k, :)');
            end
            
            H = sum1/sum2;
        end
        
        function Q = calculateQ()
            
        end
        
    end
    
end