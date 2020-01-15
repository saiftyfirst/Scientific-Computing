%=============================================================
%This function explicitly solves an ODE using the Runge-Kutta 4 method.
%
%  Input arguements:
%       f_y   = function handle of the ODE
%       y0    = initial condition 
%       dt    = time step
%       t_end = end time
%
%
%-------------------------------------------------------------

function y = rungeKatta4 (f_y, y0, dt, t_end)
    
    t = (0:dt:t_end)';
        
    looplength = length(t);
    y = zeros (1, looplength);
    y(1) = y0;
    
    for n = 1 : looplength - 1
        
        Y1     = f_y(y(n));
        Y2     = f_y(y(n) + (dt/2) * Y1);
        Y3     = f_y(y(n) + (dt/2) * Y2);
        Y4     = f_y(y(n) + dt * Y3);
        y(n+1) = y(n) + dt * (1/6) * ( Y1 + 2*Y2 + 2*Y3 + Y4 );
        
    end
    
end
