using Base.Test
import Roots.fzero

## tests:http://people.sc.fsu.edu/~jburkardt/cpp_src/test_zero/test_zero.html
function newton_baffler(x) 
    if ( x - 0.0 ) < -0.25 
        0.75 * ( x - 0 ) - 0.3125 
    elseif  ( x - 0 ) < 0.25 
        2.0 * ( x - 0 ) 
    else
        0.75 * ( x - 0 ) + 0.3125
    end
end


pathological = [
                (x ->sin ( x ) - x / 2, .1),
                (x-> 2 * x - exp ( - x ), 1),
                (x -> x * exp ( - x ), 1),
                (x -> exp ( x ) - 1 / ( 10 * x )^2, 1),
                (x -> ( x + 3 ) * ( x - 1 )^2, 1),
                
                (x -> exp ( x ) - 2 - 1 / ( 10 * x )^2 + 2 / ( 100 * x )^3, 1),
                (x -> cbrt(x), 1),
                (x -> cos ( x ) - x, 1),
                (x-> newton_baffler(x), 8),
                (x -> 20.0 * x / ( 100.0 * x^2 + 1.0), 0.1), 
                
                (x ->  ( 4.0 + x^2) * ( 2.0 + x ) * ( 2.0 - x )  / ( 16.0 * x * x * x * x + 0.00001 ), 1),
                (x -> (x == 1.0) ? 0 : sign(x-1.0) * exp(log(1e4) + log(abs(x - 1.0)) - 1.0/(x-1.0)^2), 1),
                (x -> 0.00000000001 * (x - 100.0), 1),
                (x -> 1.0 / ( ( x - 0.3 ) * ( x - 0.3 ) + 0.01 ) + 1.0 / ( ( x - 0.9 ) * ( x - 0.9 ) + 0.04 ) + 2.0 * x - 5.2, -1),
                (x -> ( 1 - 6x^2) * cbrt(x) * exp (-x^2) / (3*x), -0.25), 
                
                (x -> ( pi * ( x - 5.0 ) / 180.0 ) - 0.8 * sin ( pi * x / 180.0 ), 1),
                (x -> x^3 - 2*x - 5, 2),
                (x -> 1e6 * (x^7 -7x^6 +21x^5 -35x^4 +35x^3-21x^2+7x-1),  0.990),
                (x -> cos(100*x)-4*erf(30*x-10), 4) 
                ]
                
for (f1, x0) in pathological
    x = fzero(f1, x0)
    @assert (f1(x) == 0.0 || f1(prevfloat(x)) * f1(nextfloat(x)) <= 0 || abs(f1(x)) <= eps(float(x0))^(1/2))
end

## make a graphic comparing values
function make_graphic()
    orders = [0,1,2,5,8,16]
    N = length(orders)
    k,n = length(pathological), 50
    m = ["" for i=1:(N + 2)*k, j=1:n]
    ## how to mark values. If error we also use "x"
    function howgood(f1, xstar)
        if f1(xstar) == 0       # perfect
            "."
        elseif f1(prevfloat(xstar)) * f1(nextfloat(xstar)) < 0 # perfect
        "-"
        elseif -1e-16 < f1(xstar) < 1e-16 # approximate
            "~"
        elseif abs(xstar) > 5   # ran away yet converged
            "!"
        else                    # not close?
            "*"
        end
    end

    for i in 1:k
        f1,x = pathological[i]
        xs = linspace(-5, 5, n)
        for k in 1:n
            for j in 1:length(orders)
                try
                    xstar = fzero(f1, xs[k], order=orders[j])
                    m[(N+1)*(i-1)+j, k] = howgood(f1, xstar)
                catch e
                    m[(N+1)*(i-1)+j, k] = "x"
                end
            end
        end
    end
    println("\n----")
    print(join(mapslices(join, m, 2),"\n"))
end
