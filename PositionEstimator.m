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
        
        function W = calculateW(labels, A)
            c1 = 0;
            c2 = 0;
            for k = 2:size(labels, 2)
                c1 = c1 + (labels(:,k)*labels(:,k)');
                c2 = c2 + (labels(:,k-1)*labels(:,k)');
            end
            W = (1/size(labels, 2)-1)*(c1 - A*c2);
        end
        
        function H = calculateH(z, x, M)
            sum1 = zeros(size(z, 2), size(x, 2));
            sum2 = zeros(size(z, 2), size(x, 2));
            
            for k = 1:M
                sum1 = sum1 + (z(:, k)*x(:, k)');
                sum2 = sum2 + (x(:, k)*x(:, k)');
            end
            
            H = sum1/sum2;
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