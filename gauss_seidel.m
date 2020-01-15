function result = gauss_seidel(b, Nx, Ny)
 
    tolerence = 1e-4;
    residual = Inf;

    result = zeros((Nx + 2) * (Ny + 2), 1);
    hx_invsquared = (Nx + 1) ^ 2;
    hy_invsquared = (Ny + 1) ^ 2;
    invsq_denominator = - 2 * (hx_invsquared + hy_invsquared);
    
    Nx_padded = Nx + 2;
    Ny_padded = Ny + 2;
    
    b_compare = b(1:(ceil(Ny_padded / 2) * Nx_padded));
    
    while 1
        for j = 2:(Ny + 1)
            for k = 2:(Nx + 1)
                center_idx = (Nx_padded) * (j - 1) + k;
                result(center_idx) = ( ...
                    b(center_idx) - ...
                    ((result(center_idx - 1) + result(center_idx + 1)) * hx_invsquared) - ...
                    ((result(center_idx - Nx_padded) + result(center_idx + Nx_padded)) * hy_invsquared) ...
                ) /  invsq_denominator;                                    
            end    
        end
        
        if residual > tolerence        
            b_aux = b_compare;
            for j = 2:ceil(Ny_padded / 2)
                for k = 2:ceil(Nx_padded / 2)
                    center_idx = Nx_padded * (j - 1) + k;
                    b_aux(center_idx) = (result(center_idx) * invsq_denominator) + ...
                                        ((result(center_idx - 1) + result(center_idx + 1)) * hx_invsquared) + ...
                                        ((result(center_idx - Nx - 2) + result(center_idx + Nx + 2)) * hy_invsquared);                                 
                end    
            end
            residual = sqrt(sum(((b_compare - b_aux) .^ 2)) / length(b_compare));
        else 
            b_aux = zeros(length(b), 1);
            for j = 2:(Ny + 1)
                for k = 2:(Nx + 1)
                    center_idx = Nx_padded * (j - 1) + k;
                    b_aux(center_idx) = (result(center_idx) * invsq_denominator) + ...
                                        ((result(center_idx - 1) + result(center_idx + 1)) * hx_invsquared) + ...
                                        ((result(center_idx - Nx - 2) + result(center_idx + Nx + 2)) * hy_invsquared);                                 
                end    
            end
            residual = sqrt(sum(((b - b_aux) .^ 2)) / length(b));
            if (residual < tolerence)
                break;
            end    
        end     
    end

end