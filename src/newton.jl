# gradient based root finding methods

# convergence tests
function check_tolerance(tolerance)
    if tolerance < 0.0
        error("tolerance must be >= 0.0")
    end
end

function check_residual(fx, tolerance)
    check_tolerance(tolerance)
    abs(fx) < tolerance
end

function check_delta(delta, tolerance)
    check_tolerance(tolerance)
    abs(delta) < tolerance
end     

# Newton-Raphson method (quadratic convergence)
function newton(f::Function, fp::Function, x;
                delta_tolerance=eps(1.0),
                residual_tolerance=eps(1.0),
                max_iter=100)
    for i=1:max_iter
        fx = f(x)
        if check_residual(fx, residual_tolerance)
            return x
        end
        fpx = fp(x)
        if fpx == 0
            error("derivative is zero")
        end
        delta = fx/fpx
        x -= delta
        if check_delta(delta, delta_tolerance)
            return x
        end
    end
    error("max_iter reached without convergence")
end

# Halley's method (cubic convergence)
function halley(f::Function, fp::Function, fpp::Function, x;
                delta_tolerance=eps(1.0),
                residual_tolerance=eps(1.0),
                max_iter=100)
    for i=1:max_iter
        fx = f(x)
        if check_residual(fx, residual_tolerance)
            return x
        end
        fpx = fp(x)
        if fpx == 0
            error("derivative is zero")
        end
        fppx = fpp(x)
        if fppx == 0 # fall back to Newton
            delta = fx/fpx
        else
            delta = 2fx*fpx/(2fpx^2 - fx*fppx)
        end
        x -= delta
        if check_delta(delta, delta_tolerance)
            return x
        end
    end
    error("max_iter reached without convergence")
end

