%=============================================================
% This function evaluates the logistic growth equation for the
% given parameters.
%
%  Input arguements:
%       alpha    = equation parameter
%       beta     = equation parameter
%       p_0      = initial condition
%       t        = time steps vector        
%
%-------------------------------------------------------------

function p = logistic_growth(alpha, beta, p_0, t)
    p = beta ./ (1 + ((beta / p_0 - 1) * exp(-alpha * t)));
end
