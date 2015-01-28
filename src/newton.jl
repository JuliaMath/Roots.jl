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

function check_residual(fx, ftol)
    check_tolerance(ftol)
    abs(fx) < ftol
end

function check_delta(x, delta, xtol, xtolrel)
    check_tolerance(xtol)
    abs(delta) <= max(xtol, abs(x)*xtolrel)
end     


## Order 1 secant method
function secant_method(f::Function, x0::Real, x1::Real;
                       ftol::Real     = 10.0 * eps(one(eltype(float(x1)))),
                       xtol::Real     = 10.0 * eps(one(eltype(float(x1)))),
                       xtolrel::Real = eps(one(eltype(float(x1)))),
                       maxeval::Int   = 100, 
                       verbose::Bool  = false,
                kwargs...)

    a,b = sort([x0,x1])
    fa, fb = f(a), f(b)

    abs(fb) <= ftol && return(b)
    abs(b-a) <= max(xtol, abs(b)*xtolrel) && throw(ConvergenceFailed("a,b too close"))
    
    try
        abs(fb) <= ftol && throw(StateConverged(b))

        for i in 1:maxeval
            verbose && println("$i: a=$a, b=$b, f(b)=$fb")
            sec = secant(fa, fb, a, b)
            ((sec == 0.0) | isnan(sec)) && throw(ConvergenceFailed("Division by 0"))
            inc = fb / sec
            
            a, b = b, b-inc
            fa, fb = fb, f(b)


            (isinf(b) | isinf(fb)) && throw(ConvergenceFailed("Function values diverged"))
            abs(fb) <= ftol && throw(StateConverged(b))
            
        end
        throw(ConvergenceFailed("More than $maxeval steps taken before convergence"))
        
    catch e
        if isa(e, StateConverged)
            e.x0
        else
            rethrow(e)
        end
    end
end


# Newton-Raphson method (quadratic convergence)
function newton(f::Function, fp::Function, x;
                ftol::Real       = 100 * eps(one(eltype(float(x)))),
                xtol::Real       = 10  * eps(one(eltype(float(x)))),
                xtolrel::Real    = eps(one(eltype(float(x)))),
                maxeval::Integer = 100,
                verbose::Bool    = false
                )

    cvg = false
    for i=1:maxeval
        fx = f(x)
        if check_residual(fx, ftol)
            cvg = true
            break
        end

        fpx = fp(x)
        if fpx == 0
            throw("derivative is zero")
        end

        del = fx/fpx
        x -= del

        if check_delta(x, del, xtol, xtolrel)
            cvg = true
            break
        end
        verbose && println("xn = $x, f(xn) = $(f(x)), step=$i")
    end
    cvg || error("$maxeval steps taken without convergence")
    
    return x
end



# Halley's method (cubic convergence)
function halley(f::Function, fp::Function, fpp::Function, x;
                ftol::Real       = 100 * eps(one(eltype(float(x)))),
                xtol::Real       = 10  * eps(one(eltype(float(x)))),
                xtolrel::Real    = eps(one(eltype(float(x)))),
                maxeval::Integer = 100,
                verbose::Bool    = false
)

    cvg = false
    for i=1:maxeval
        fx = f(x)
        if check_residual(fx, ftol)
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
        if check_delta(x, del, xtol, xtolrel)
            cvg = true
        end
        verbose && println("xn = $x, f(xn) = $(f(x)), step=$i")
    end
    cvg || throw("$maxeval steps taken without convergence")
    
    return x

end

