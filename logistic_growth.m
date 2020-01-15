function p = logistic_growth(alpha, beta, p_0, t)
    p = beta ./ (1 + ((beta / p_0 - 1) * exp(-alpha * t)));
end