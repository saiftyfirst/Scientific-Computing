function b_padded = vector_padding(b, Nx, Ny)
    b_padded = zeros((Nx + 2) * (Ny + 2), 1);
    for j = 2:(Ny + 1)
        b_idx = ((j - 2) * Nx) + 1;
        b_padded_idx = (Nx + 2) * (j - 1) + 2;
        b_padded(b_padded_idx:(b_padded_idx + Nx - 1)) = b(b_idx:(b_idx + Nx - 1));
    end
end