format short
% PLOTTING ANALYTIC SOLUTION OF:
%  p(t) = 10 / (1 + 9 * e ^ (-t)) 
% Verhulst logistic growth with parameters:
beta = 10; alpha = 1;
p_0 = 1;                %initial condition
time = linspace(0,20);  %time vector
population = logistic_growth(alpha, beta, p_0, time); %evaluation of function
plot(time, population);                               %plot function
title("Plot of p(t) = 10 / (1 + 9 * e ^ (-t)))");     %plot title
axis([0 12 0 18]);                                    %adjusting plot axes
xlabel('t')                                           %label to X axis
ylabel('y(t)')                                        %label to Y axis


% NUMERICAL METHODS
%Given Data / Instructions
p_dot = @(p)(1 - p / 10) * p;                         %function handle for ODE  
p_analytical = @(t) 10 ./ ( 1 + 9 * exp(-t));         %function handle for analytic solution

p_0 = 1;                    %initial condition
dt = [1/8 1/4 1/2 1];       %time steps to solve for

t_end = 5;                  %end time
method_labels = ["Analytical", "Euler", "Heun", "RK4"];

% Outer loop: Loops through each method used
% Inner loop: Call appropriate method and plots it against the exact
% solution for each variation of dt
for i = 1:length(method_labels) - 1
    figure('Name', method_labels(i+1));
    rel_error_vector = zeros(1, length(dt));
    approx_error_vector = zeros(1, length(dt));
    
    for j = 1:length(dt)
        % evaluation
        time = 0:dt(j):t_end;
        population_approx = approximate_diff(p_dot, p_0, dt(j), t_end, i);
        population_exact = p_analytical(time);        
        
        % error calculation 
        rel_error_vector(j) = errorRelative(population_exact, ...
                              population_approx, dt(j), t_end);
        if j ~= 1
            approx_error_vector(j) = errorApprox(best_approx, ... 
                                     population_approx, ...
                                     dt(1), dt(j), t_end);
        else 
            best_approx = population_approx;
        end
                         
        % plotting
        subplot(2, 2, j);  
        hold on;
        plot(time, population_exact);
        plot(time, population_approx);
        xlabel('t'); ylabel('p(t)')
        axis([0 (t_end + 0.5) 0 15]);
        title(strcat('dt=', string(dt(j))));
        legend(method_labels(1), method_labels(i + 1));
    
    end
    
    print_table(method_labels(i + 1), flip(dt), flip(rel_error_vector), flip(approx_error_vector));
    
end

function print_table(name, dt, rel_error_vector, approx_error_vector) 
    disp(name)
    colNames = strcat("dt=",string(dt));
    rowNames = {'error', 'error red.', 'error app.'};
    
    data = [rel_error_vector; ...
            [0 rel_error_vector(2:end)] ./ [1 rel_error_vector(1:end-1)]; ... 
            approx_error_vector];
    array2table(data, 'RowNames', rowNames, 'VariableNames', colNames)
end

%Function to switch through the methods in the above loop
function y = approximate_diff(f_y, y0, dt, t_end, mode) 
    switch mode
        case 1
            y = explicitEuler(f_y, y0, dt, t_end);
        case 2
            y = heunMethod(f_y, y0, dt, t_end);
        case 3
            y = rungeKatta4(f_y, y0, dt, t_end);
        otherwise
            warning('Stop choosing non-existing methods!');
    end
end



