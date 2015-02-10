## A collection of historical methods for pedagogical purposes.
##
## secant_method, newton, and halley
##
## These have an argument `verbose` that can be specified to get a trace of the algorithm



## Order 1 secant method
function secant_method(f::Function, x0::Real, x1::Real;
                       ftol::Real     = 10.0 * eps(one(eltype(float(x1)))),
                       xtol::Real     = 10.0 * eps(one(eltype(float(x1)))),
                       xtolrel::Real = eps(one(eltype(float(x1)))),
                       maxeval::Int   = 100, 
                       verbose::Bool  = false,
                kwargs...)

    if (ftol < 0) | (xtol < 0) | (xtolrel < 0)
        throw(ConvergenceFailed("Tolerances must be non-negative"))
    end
    
    fx0, fx1 = f(x0), f(x1)

    abs(fx1) <= ftol && return(x1)
    abs(x1 - x0) <= max(xtol, abs(x1)*xtolrel) && throw(ConvergenceFailed("x0, x1 too close"))
    
    try
        abs(fx1) <= ftol && throw(StateConverged(x1))

        for i in 1:maxeval
            verbose && println("$i: x_$(i-1)=$x0, x_$i=$x1, f(x_$i)=$fx1")
            appsecinv = (x1 - x0) / (fx1 - fx0)
            ((appsecinv == 0.0) | isnan(appsecinv)) && throw(ConvergenceFailed("Division by 0"))

            inc = fx1 * appsecinv
            x0, x1 = x1, x1 - inc
            fx0, fx1 = f(x0), f(x1)


            abs(fx1) <= ftol && throw(StateConverged(x1))
            (isinf(x1) | isinf(fx1)) && throw(ConvergenceFailed("Function values diverged"))
            
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

    if (ftol < 0) | (xtol < 0) | (xtolrel < 0)
        throw(ConvergenceFailed("Tolerances must be non-negative"))
    end
    
    cvg = false
    for i=1:maxeval
        fx = f(x)
        verbose && println("x_$(i-1) = $x, f(x_$(i-1)) = $(f(x))")
        
        if abs(fx) <= ftol
            cvg = true
            break
        end

        fpx = fp(x)
        if fpx == 0
            throw("derivative is zero")
        end

        del = fx/fpx
        x -= del

        if abs(del) <= max(xtol, abs(x)*xtolrel)
            cvg = true
            break
        end
        

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

    if (ftol < 0) | (xtol < 0) | (xtolrel < 0)
        throw(ConvergenceFailed("Tolerances must be non-negative"))
    end

    cvg = false
    for i=1:maxeval
        fx = f(x)
        verbose && println("x_$(i-1) = $x, f(x_$(i-1)) = $(f(x))")
        
        if abs(fx) <= ftol
            cvg = true
            break
        end
        fpx = fp(x)
        fpx == 0 && throw(ConvergenceFailed("Derivative is zero"))

        fppx = fpp(x)
        if fppx == 0 # fall back to Newton
            del = fx/fpx
        else
            del = 2fx*fpx/(2fpx^2 - fx*fppx)
        end
        x -= del

        if abs(del) <= max(xtol, abs(x) * xtolrel)
            cvg = true
        end
    end
    cvg || throw("$maxeval steps taken without convergence")
    
    return x

end

