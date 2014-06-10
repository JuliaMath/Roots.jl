## A collection of historical methods for pedagogical purposes.
##
## secant_method, newton, and halley
##
## These have an argument `verbose` that can be specified to get a trace of the algorithm

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


## Order 1 secant method
function secant_method(f::Function, x0::Real, x1::Real;
                tol::Real   = 10.0 * eps(one(eltype(float(x1)))),
                delta::Real =  zero(x1),
                max_iter::Int=100, 
                verbose::Bool=false,
                kwargs...)

    a, b, fa, fb = x0, x1, f(x0), f(x1)
    
    try
        fb == 0 && throw(StateConverged(b))

        for i in 1:max_iter
            verbose && println("a=$a, b=$b, ctr=$(i-1)")

            inc = fb/secant(fa, fb, a, b)
            a, b = b, b-inc
            fa, fb = fb, f(b)
            abs(inc) <= delta && throw(StateConverged(b))
            abs(fb) <= tol && throw(StateConverged(b))

        end

        throw(ConvergenceFailed())

    catch e
        if isa(e, StateConverged)
            e.x0
        else
            throw(e)
        end
    end
end


# Newton-Raphson method (quadratic convergence)
function newton(f::Function, fp::Function, x;
                delta::Real       = 10  * eps(one(eltype(float(x)))),
                tol::Real         = 100 * eps(one(eltype(float(x)))),
                max_iter::Integer = 100,
                verbose::Bool     = false
                )

    cvg = false
    for i=1:max_iter
        fx = f(x)
        if check_residual(fx, tol)
            cvg = true
            break
        end

        fpx = fp(x)
        if fpx == 0
            throw("derivative is zero")
        end

        del = fx/fpx
        x -= del

        if check_delta(del, delta)
            cvg = true
            break
        end
        verbose && println("xn = $x, f(xn) = $(f(x)), step=$i")
    end
    cvg || error("$max_iter steps taken without convergence")
    
    return x
end



# Halley's method (cubic convergence)
function halley(f::Function, fp::Function, fpp::Function, x;
                delta::Real       = 10  * eps(one(eltype(float(x)))),
                tol::Real         = 100 * eps(one(eltype(float(x)))),
                max_iter::Integer = 100,
                verbose::Bool     = false
)

    cvg = false
    for i=1:max_iter
        fx = f(x)
        if check_residual(fx, tol)
            cvg = true
            break
        end
        fpx = fp(x)
        if fpx == 0
            throw("derivative is zero")
        end
        fppx = fpp(x)
        if fppx == 0 # fall back to Newton
            del = fx/fpx
        else
            del = 2fx*fpx/(2fpx^2 - fx*fppx)
        end
        x -= del
        if check_delta(del, delta)
            cvg = true
        end
        verbose && println("xn = $x, f(xn) = $(f(x)), step=$i")
    end
    cvg || throw("$max_iter steps taken without convergence")
    
    return x

end

