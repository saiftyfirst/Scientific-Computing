%=============================================================
%This function explicitly solves an ODE using the Euler method.
%
%  Input arguements:
%       f_y   = function handle of the ODE
%       y0    = initial condition 
%       dt    = time step
%       t_end = end time
%
%       y(n+1)    = y(n) + dt * (f(y(n)))
%       
%
%-------------------------------------------------------------

function y = explicitEuler (f_y, y0, dt, t_end)
    
    t = (0:dt:t_end)';
    fy = f_y(y0); %evaluate f_y at initial y0
    
    looplength = length(t);
    y = zeros(1, looplength);
    y(1) = y0;
    
    for n=1:looplength-1
        y(n+1) = y(n)+dt*fy;
        fy = f_y(y(n+1));
    end
    
    
end
