format short
% % PLOTTING ANALYTIC SOLUTION OF:
% %  p(t) = 200 / (20 - 1- * e ^ (-7*t))
% % Verhulst logistic growth with parameters:
% beta = 10; alpha = 7;
% p_0 = 20;                %initial condition
% time = linspace(0,8);  %time vector
% population = logistic_growth(alpha, beta, p_0, time); %evaluation of function
% plot(time, population);                               %plot function
% title("Plot of p(t) = 200 / (20 - 1- * e ^ (-7*t))");     %plot title
% axis([0 8 0 24]);                                    %adjusting plot axes
% xlabel('t')                                           %label to X axis
% ylabel('p(t)')      

% common parameters for all functions
p_dot = @(p) 7 * (1 - p / 10) * p;  % function handle for ODE  
p_analytical = @(t) 200 ./ ( 20 - 10 * exp(-7 * t)); %function handle for analytic solution
p_0 = 20;
dt = [1/32 1/16 1/8 1/4 1/2];
legend = ['Analytic' strcat('dt=', string(dt))];
t_end = 5;

do_work(p_dot, p_analytical, p_0, dt, t_end, @approximate_diff, 1, 'Explicit Euler', legend);
do_work(p_dot, p_analytical, p_0, dt, t_end, @approximate_diff, 2, 'Huen', legend);
do_work(p_dot, p_analytical, p_0, dt, t_end, @approximate_diff, 4, 'Implicit Euler', legend);
do_work(p_dot, p_analytical, p_0, dt, t_end, @approximate_diff, 5, '2nd Order Adams Moulton', legend);
do_work(p_dot, p_analytical, p_0, dt, t_end, @approximate_diff_wk3, 1, '1st Linear Adams Moulton', legend);
do_work(p_dot, p_analytical, p_0, dt, t_end, @approximate_diff_wk3, 2, '2nd Linear Adams Moulton', legend);

function do_work(p_dot, p_analytical, p_0, dt, t_end, function_group, method, title_, legend_)
    figure('Name', title_);
    hold on;
    time = 0:dt(1):t_end;
    population_exact = p_analytical(time);        
    plot(time, population_exact);
    approx_error_vector = zeros(1, length(dt));
    for j = 1:length(dt)
        time = 0:dt(j):t_end;
        % exact solution
        population_exact = p_analytical(time);        
        % approximate 
        population_approx = function_group(p_dot, p_0, dt(j), t_end, method);
        % clean approximation data
        population_approx(population_approx == -inf) = min(population_approx(isfinite(population_approx)));
        population_approx(population_approx == inf) = max(population_approx(isfinite(population_approx)));
        % slice time in case of newton method failure
        time = time(1:length(population_approx));
        approx_error_vector(j) = error_approx(population_exact, population_approx, dt(j), t_end);
        plot(time, population_approx);
    end
    xlabel('t'); ylabel('p(t)');
    axis([0 5 0 20]);
    title(title_);
    legend(legend_);
    print_table(title_, flip(dt), flip(approx_error_vector));
end

function print_table(name, dt, approx_error_vector) 
    disp(name)
    colNames = strcat("dt=",string(dt));
    rowNames = {'error', 'error red.'};
    error_red = [0 approx_error_vector(1:end-1)] ./ [1 approx_error_vector(2:end)];
    error_red(isnan(error_red)) = 0;
    data = [approx_error_vector; error_red];
    array2table(data, 'RowNames', rowNames, 'VariableNames', colNames)
end

