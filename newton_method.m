function x_approx = newton_method(G, x_0) 
    TOLERENCE = 0.0001;
    MAX_ITERATIONS = 20;
    eval = G(x_0);
    x_approx = x_0;
    iterations = 0;
    while abs(eval) > TOLERENCE 
        x_approx = x_approx - (eval / approximate_slope(G, x_approx));
        eval = G(x_approx);
        iterations = iterations + 1;
        if iterations > MAX_ITERATIONS
            x_approx = nan;
            return;
        end
    end
end

function m = approximate_slope(G, x)
    SLOPE_APPROX_DX = 0.001;
    m = (G(x + SLOPE_APPROX_DX) - G(x)) / SLOPE_APPROX_DX;
end
