N = [3 7 15 31 63 127];

[full_direct_time, sparse_direct_time, g_seidel_time, ...
    full_direct_storage, sparse_direct_storage, g_seidel_storage, errors] = ...
    deal(zeros(1, length(N)), zeros(1, length(N)), zeros(1, length(N)), ...
        zeros(1, length(N)), zeros(1, length(N)), zeros(1, length(N)), zeros(1, length(N)));

analytic_solution = @(x, y) (sin(pi * x)) .* (sin(pi * y));

for i = 1:length(N)
        
    fig = figure('NumberTitle', 'off', 'Name', sprintf("N=" + string(N(i))));
    
    [A, b] = discretize(N(i), N(i));
    
    range = 0:(1 / (N(i) + 1)):1;
    [X, Y] = meshgrid(range, range);
    
    tic;
    solution = A \ b;
    full_direct_time(i) = toc;
    full_direct_storage(i) = numel(A) + length(b);
    fd_padded_solution = vector_padding(solution, N(i), N(i));
    fd_padded_solution = transpose(reshape(fd_padded_solution, [(N(i) + 2), (N(i) + 2)]));
    % plotting direct solutions
    subplot(3,2,1);
    surf(X, Y, fd_padded_solution);
    title('Full Direct');  
    colorbar;
    subplot(3,2,2);
    contour(X, Y, fd_padded_solution);
    title('Full Direct');  
    colorbar;
    
    A_sparse = sparse(A);
    tic;
    solution = A_sparse \ b;
    sparse_direct_time(i) = toc;
    sparse_direct_storage(i) = length(A) * 5 + length(b);
    sd_padded_solution = vector_padding(solution, N(i), N(i));
    sd_padded_solution = transpose(reshape(sd_padded_solution, [(N(i) + 2), (N(i) + 2)]));    
    % plotting sparse direct solutions
    subplot(3,2,3);
    surf(X, Y, sd_padded_solution);
    title('Sparse Direct');
    colorbar;
    subplot(3,2,4);
    contour(X, Y, sd_padded_solution);    
    title('Sparse Direct');
    colorbar;
    
    b_padded = vector_padding(b, N(i), N(i));
    tic;
    solution = gauss_seidel(b_padded, N(i), N(i));
    g_seidel_time(i) = toc;
    shaped_solution = transpose(reshape(solution, [(N(i) + 2), (N(i) + 2)]));
    g_seidel_storage(i) = length(b_padded) * 2;
    
    % plotting gauss-seidel solutions
    subplot(3,2,5);
    surf(X, Y, shaped_solution);
    title('Gauss-Seidel');        
    colorbar;   
    subplot(3,2,6);
    contour(X, Y, shaped_solution);
    title('Gauss-Seidel');        
    colorbar;   
    
    analytic_results = analytic_solution(X(:), Y(:));
    errors(i) = sqrt((sum((analytic_results - solution) .^ 2) / (N(i) ^ 2)));

end    

print_table("Full Matrix Direct Solver", N, [full_direct_time; full_direct_storage]);
print_table("Sparse Matrix Direct Solver", N, [sparse_direct_time; sparse_direct_storage]);
print_table("Gauss-Seidel Solver", N, [g_seidel_time; g_seidel_storage]);
print_error_table(N, errors);

function print_table(name, N, data) 
    disp(name)
    colNames = strcat("Nx,Ny = ", string(N));
    rowNames = {'runtime', 'storage'};
    array2table(data, 'RowNames', rowNames, 'VariableNames', colNames)
end

function print_error_table(N, error_data) 
    disp("Error Analysis")
    colNames = strcat("Nx,Ny = ", string(N));
    rowNames = {'error', 'error_reduction'};
    reduction_data = [0 error_data(1:end-1)] ./ [1 error_data(2:end)];
    data = [error_data; reduction_data];
    array2table(data, 'RowNames', rowNames, 'VariableNames', colNames)
end
