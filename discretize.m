function [A, b] = discretize(Nx, Ny)
    dimensions_product = Nx * Ny;
    
    hx = 1 / (Nx +1);
    hy = 1 / (Ny +1);
    
    hx_invsqr = (1 / hx) ^ 2;
    hy_invsqr = (1 / hy) ^ 2;
    
    % off-diagonals
    one_translated_diagonal = setdiff(1:(dimensions_product - 1), Nx:Nx:(dimensions_product - 1));
    
    row_vector = [1:(dimensions_product - Nx) one_translated_diagonal];
    col_vector = [(Nx + 1):dimensions_product (one_translated_diagonal + 1)];
    value_vector = [(zeros(1, dimensions_product - Nx) + hy_invsqr) (zeros(1, length(one_translated_diagonal)) + hx_invsqr)];

    A = sparse(row_vector, col_vector, value_vector, dimensions_product, dimensions_product);
    A = A + transpose(A);
    
    % diagonals
    row_vector = 1:dimensions_product;
    col_vector = 1:dimensions_product;
    value_vector = zeros(1, dimensions_product) + (-2 * (hx_invsqr + hy_invsqr));
    
    A = A + sparse(row_vector, col_vector, value_vector, dimensions_product, dimensions_product);

%     
%     % off diagonals (Nx translated)
%     row_vector = [row_vector 1:(dimensions_product - Nx) (Nx + 1):dimensions_product];
%     col_vector = [col_vector (Nx + 1):dimensions_product 1:(dimensions_product - Nx)];
%     value_vector = [value_vector (zeros(1, (dimensions_product - Nx) * 2) + hy_invsqr)];
%     
%     % off diagonals (1 translated)
%     row_vector = [row_vector 1:(dimensions_product - 1)];
%     col_seq = setdiff(1:(dimensions_product - 1), Nx:Nx:(dimensions_product - 1)) + 1;
%     col_vector = [col_vect col_sequence];
%     
%     row_vector = [(1:dimensions_product - 1) 1:(dimensions_product - Nx)];
%     col_vector = [2:(dimensions_product) (Nx + 1):dimensions_product];
%     value_vector = [(zeros(1, (dimensions_product - 1)) + hx_invsqr) (zeros(1, (dimensions_product - Nx)) + hy_invsqr)];
%     A = sparse(row_vector, col_vector, value_vector, dimensions_product, dimensions_product);
%     A = A + transpose(A);
%     

    
    b = ones(dimensions_product, 1);
  
end