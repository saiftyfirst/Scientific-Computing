%=============================================================
% This function calculates the error between the the best 
% numerical solution and another numeric solution.
%
%  Input arguements:
%       yBest     = best numeric solution vector
%       yTest     = numerical solution vector we want to test
%       dt        = time step of test solution
%       dt_best   = time step of best solution        
%       t_end     = end time
%
%       error    = .......
%
%-------------------------------------------------------------

function err = error_relative (y_best, y_test, dt_best, dt, t_end)
    
    skip = dt / dt_best;                               %calculates how many values of the yBest vector to skip to create the vector comparable to the yTest
    yBest = y_best(1:skip:length(y_best));               %creates the yBest vector with only the values that will be compared with yTest
    
    err = sqrt((dt / t_end) * sum((y_test - yBest) .^ 2)); %calculates the relative error with the given formula, in a vectorized format
    
end


