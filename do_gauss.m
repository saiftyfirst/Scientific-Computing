% Iteration 1: Only square, non-singular matrices.
function x = do_gauss(A, b)
    n = length(A);
    [At, c] = forward_eliminate(A,b, n);
    x = backward_substitute(At, c, n);
end

function [At, c] = forward_eliminate(A, b, n) 

    At = zeros(n);
    c = zeros(n, 1);
    
    for j = 1:n-1 % for each row in A

        % Store already finalised elememts into At
        for k = j:n
            At(j,k) = A(j,k);
        end
        c(j) = b(j);
    
        % Elimination from the following column
        for i = j+1:n
            factor = A(i,j) / At(j,j);
            for k = j+1:n
                A(i,k) = A(i,k) - factor * At(j,k);
            end
            b(i) = b(i) - factor * c(j);
        end
        
    end
    
    c(n) = b(n);
    At(n,n) = A(n,n);
    
end

function x = backward_substitute(At, c, n) 
    x = zeros(n, 1);
    for i = n:-1:1
     for j = n:-1:i+1
        c(i) = c(i) - At(i,j) * x(i+1);     
     end
     x(i) = c(i) / At(i, i);
    end
end