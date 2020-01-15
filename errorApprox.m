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

function err = errorApprox (yBest, yTest, dt_best, dt, t_end)
    
    skip = dt / dt_best;                               %calculates how many values of the yBest vector to skip to create the vector comparable to the yTest
    yBest = yBest(1:skip:length(yBest));               %creates the yBest vector with only the values that will be compared with yTest
    
    err = sqrt( (dt/t_end) * sum( (yTest-yBest).^2) ); %calculates the relative error with the given formula, in a vectorized format
    
end


