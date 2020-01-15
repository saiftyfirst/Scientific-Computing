%=============================================================
%This function explicitly solves an ODE using the Heun method.
%
%  Input arguements:
%       f_y   = function handle of the ODE
%       y0    = initial condition 
%       dt    = time step
%       t_end = end time
%
%       y(n+1)    = y(n) + dt * (1/2) * (f(y(n)) + f(y(n+1)))
%       f(y(n+1)) = y(n) + dt * f(y(n)) 
%
%-------------------------------------------------------------

function y = heunMethod (f_y, y0, dt, t_end)
    
    t = (0:dt:t_end)'; %create time vector (is it needed...?)
    
    looplength = length(t); %stopping point for loop
    y = zeros(1,looplength); %initialize y vector
    fy = zeros(2,1);
    y(1) = y0;


    
    for n=1:length(t)-1
        fy(1) = f_y(y(n));
        fy(2) = f_y(y(n) + dt * fy(1));
        y(n+1)    = y(n) + dt * (1/2) * (fy(1) + fy(2));
           
    end
    
    
end
