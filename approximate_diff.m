% Solves initial value problems where y(0) = y_0
function y = approximate_diff(f, y_0, dt, t_end, mode)
    switch mode
        case 1
            y = explicit_eulerify(f, y_0, dt, t_end);
        case 2
            y = heunify(f, y_0, dt, t_end);
        case 3
            y = runge_kutta_baby(f, y_0, dt, t_end);
        case 4
            y = implicit_eulerify(f, y_0, dt, t_end);
        case 5
            y = adams_moulton_second(f, y_0, dt, t_end);
        otherwise
            warning('Stop choosing non-existing methods!');
    end
end

function y = explicit_eulerify(f, y_0, dt, t_end) 
   y = [y_0, zeros(1, t_end/dt)];
   for i = 2:length(y)
       y(i) = y(i-1) + f(y(i-1)) * dt;
   end
end

function y = heunify(f, y_0, dt, t_end) 
   y = [y_0, zeros(1, t_end/dt)];
   for i = 2:length(y)
       y(i) = y(i-1) + (f(y(i-1)) + f(y(i-1) + f(y(i-1)) * dt))* 0.5 * dt;
   end
end

function y = runge_kutta_baby(f, y_0, dt, t_end)
   y = [y_0, zeros(1, t_end/dt)];
   for i = 2:length(y)
       k = [f(y(i-1)), zeros(1, 3)];
       for j = 2:4
           k(j) = f(y(i-1) + k(j-1) * round(2 * (j / 4 - 0.1)) * dt);
       end
       y(i) = y(i-1) + (sum([1 2 2 1 .* x])) / 6 * dt;
   end
end

function y = implicit_eulerify (f, y_0, dt, t_end)
   y = [y_0, zeros(1, t_end/dt)];    
   for i = 2:length(y)
       F = @(Y) Y - f(Y) * dt - y(i-1);
       y(i) = newton_method(F, y(i-1));
       if isnan(y(i))
           y = y(1:i-1);
           return;
       end
   end
end


function y = adams_moulton_second (f, y_0, dt, t_end)
   y = [y_0, zeros(1, t_end/dt)];    
   for i = 2:length(y)
       F = @(Y) Y - (f(Y) + f(y(i-1))) * dt * 0.5 - y(i-1);
       y(i) = newton_method(F, y(i-1));
       if isnan(y(i))
           y = y(1:i-1);
           return;
       end
   end
end
