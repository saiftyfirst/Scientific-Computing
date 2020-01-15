function [A, b] = discretize(Nx, Ny)
    dimensions_product = Nx * Ny;
    
    hx = 1 / (Nx +1);
    hy = 1 / (Ny +1);
    
    hx_invsqr = (1 / hx) ^ 2;
    hy_invsqr = (1 / hy) ^ 2;
    
    x_vector = hx_invsqr + zeros(dimensions_product -1, 1);
    x_vector(Nx:Nx:length(x_vector)) = 0;
    
    main_diagonal = diag((-2 * (hx_invsqr + hy_invsqr)) + zeros(dimensions_product, 1));
    y_diagonal = diag(hy_invsqr + zeros(dimensions_product - Nx, 1), Nx);
    x_diagonal = diag(x_vector, 1);
    
    A = y_diagonal + x_diagonal;
    A = A + transpose(A) + main_diagonal;
   
    f = @(x, y) -2 * (pi ^ 2) * sin(pi * x) .* sin(pi * y);
    
    range_horizontal = repelem(hx:hx:1-hx, length(hy:hy:1-hy));
    range_vertical = repmat(hy:hy:1-hy, 1, length(hx:hx:1-hx));
    
    b = transpose(f(range_horizontal, range_vertical));
  
end