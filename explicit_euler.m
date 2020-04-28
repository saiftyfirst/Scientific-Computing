function T = explicit_euler(Nx, Ny, dt, t_end)
    T = ones(Nx * Ny, 1);
    T = vector_padding(T, Nx, Ny);

    hx_invsquared = (Nx + 1) ^ 2;
    hy_invsquared = (Ny + 1) ^ 2;
    invsq_denominator = -2 * (hx_invsquared + hy_invsquared);
    
    for i = 0:dt:t_end
        T = iterative_update(T, dt, Nx, Ny, hx_invsquared, hy_invsquared, invsq_denominator);
    end    
end

function T = iterative_update(T_prev, dt, Nx, Ny, hx_invsquared, hy_invsquared, invsq_denominator)
    T = T_prev;
    
    Nx_padded = Nx + 2;
    
    for j = 2:(Ny + 1)
        for k = 2:(Nx + 1)
            center_idx = (Nx_padded) * (j - 1) + k;
            T(center_idx) = T_prev(center_idx) + dt * ((T_prev(center_idx) * invsq_denominator) + ...
                ((T_prev(center_idx - 1) + T_prev(center_idx + 1)) * hx_invsquared) + ...
                ((T_prev(center_idx - Nx_padded) + T_prev(center_idx + Nx_padded)) * hy_invsquared));                                    
        end    
    end
end
