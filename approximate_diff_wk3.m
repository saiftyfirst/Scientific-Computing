function y = approximate_diff_wk3(~, y_0, dt, t_end, mode)
    switch mode
        case 1
            y = adams_moulton_linearise_uno(y_0, dt, t_end);
        case 2
            y = adams_moulton_linearise_duo(y_0, dt, t_end);
        otherwise
            warning('Stop choosing non-existing methods!');
    end
end

function y = adams_moulton_linearise_uno(y_0, dt, t_end) 
   f = @(y_prev, dt) (y_prev * (1 + 7 * dt) - (7 * y_prev * y_prev * dt / 20)) / ... 
                        (1 + 7 * y_prev * dt / 20);
   y = [y_0, zeros(1, t_end/dt)];
   for i = 2:length(y)
       y(i) = f(y(i-1), dt);
   end
end

function y = adams_moulton_linearise_duo(y_0, dt, t_end) 
   f = @(y_prev, dt) (y_prev + dt * 0.5 * 7 * y_prev * (1 - y_prev / 10)) / ...
                        (1 - dt * 0.5 * 7 * y_prev * (1 - y_prev / 10));
   y = [y_0, zeros(1, t_end/dt)];
   for i = 2:length(y)
       y(i) = f(y(i-1), dt);
   end
end