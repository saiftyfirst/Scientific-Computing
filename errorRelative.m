%=============================================================
% This function calculates the error between the exact (analytic)
% and the approximate (numerical) solution.
%
%  Input arguements:
%       y     = analytic solution vector
%       yA    = numerical solution vector
%       dt    = time step
%       t_end = end time
%
%       error    = .......
%
%-------------------------------------------------------------

function err = errorRelative (y, yN, dt, t_end)
    
    err = sqrt( (dt/t_end) * sum( (y-yN).^2) );     %calculates the relative error with the given formula, in a vectorized format
    
end

