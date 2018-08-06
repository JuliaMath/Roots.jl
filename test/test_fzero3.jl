import Roots.fzero

# test robustness of derivative free algorithms


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



mutable struct CallableFunction
f
n
end
(F::CallableFunction)(x) = (F.n+=1; F.f(x))

pathological = [
                (x ->sin( x ) - x / 2, .1),
                (x-> 2 * x - exp( - x ), 1),
                (x -> x * exp( - x ), 1),
                (x -> exp( x ) - 1 / ( 10 * x )^2, 1),
                (x -> ( x + 3 ) * ( x - 1 )^2, 1),
                
                (x -> exp( x ) - 2 - 1 / ( 10 * x )^2 + 2 / ( 100 * x )^3, 1),
                (x -> cbrt(x), 1),
                (x -> cos( x ) - x, 1),
                (x-> newton_baffler(x), 8),
                (x -> 20.0 * x / ( 100.0 * x^2 + 1.0), 0.095), 
                
                (x ->  ( 4.0 + x^2) * ( 2.0 + x ) * ( 2.0 - x )  / ( 16.0 * x * x * x * x + 0.00001 ), 1),
                (x -> (x == 1.0) ? float(0) : sign(x-1.0) * exp(log(1e4) + log(abs(x - 1.0)) - 1.0/(x-1.0)^2), 1),
                (x -> 0.00000000001 * (x - 100.0), 1),
                (x -> 1.0 / ( ( x - 0.3 ) * ( x - 0.3 ) + 0.01 ) + 1.0 / ( ( x - 0.9 ) * ( x - 0.9 ) + 0.04 ) + 2.0 * x - 5.2, -1),
                (x -> ( 1 - 6x^2) * cbrt(x) * exp(-x^2) / (3*x), -0.25), 
                
                (x -> ( pi * ( x - 5.0 ) / 180.0 ) - 0.8 * sin( pi * x / 180.0 ), 1),
                (x -> x^3 - 2*x - 5, 2),
                (x -> 1e6 * (x^7 -7x^6 +21x^5 -35x^4 +35x^3-21x^2+7x-1),  0.990),
                (x -> cos(100*x)-4*erf(30*x-10), 0.0) 
                ]
      

function howgood(f1, xstar)
    if isnan(xstar)
        "x"
    elseif iszero(f1(xstar))       # perfect
        "."
    elseif f1(prevfloat(xstar)) * f1(nextfloat(xstar)) < 0 # perfect
        "."
    elseif abs(f1(xstar)) < 8 * eps(xstar) # nearly perfect
        "."
    elseif abs(f1(xstar)) <= 10eps(xstar)^(1/3)
        "~"
    elseif abs(xstar) > 5   # ran away yet converged
        "!"
    else                    # not close?
        "*"
    end
end

## make a graphic comparing values
function run_robustness_test(bigfloats=false)
    orders = [Order0(),Order1(),Order2(),Order5(),Order8(),Order16()]
    N = length(orders) + 1
    k,n = length(pathological), 50
    m = ["" for i=1:(N + 2)*k, j=1:n]
    ## how to mark values. If error we also use "x"
  
    marks = ["HP", "0 ", "1 ", "2 ", "5 ", "8 ", "16"]
    io = IOBuffer()

    for i in 1:k

        println("")
        println("-- $i ---")
        println("")


        f1,x = pathological[i]
        if bigfloats
            xs = -big(5):10/50:5
        else
            xs = linspace(-5, 5, n)
        end

        calls = Int[0]
        hits = 0
        j = 1
        print(io, marks[j], " ")
        for k in 1:n
            F = CallableFunction(f1, 0)
            xstar = Roots.dfree(F, xs[k])
            push!(calls, F.n)
            val = howgood(f1, xstar)
            val in [".", "~"] && (hits += 1)
            print(io, val)
        end
        print(String(take!(io)))
        print(" h=$hits n = $(sum(calls)) max=$(maximum(calls))")
        println("")
        
        for j in 2:N
            calls = Int[0]
            hits = 0
            print(io, marks[j], " ")
            for k in 1:n
                M = orders[j-1]
                F = CallableFunction(f1, 0)
                state = Roots.init_state(M, F, xs[k])
                options = Roots.init_options(M, state)
                tracks = Roots.NullTracks()
                xstar = find_zero(M, F, options, state, tracks)
                push!(calls, F.n)
                val = howgood(f1, xstar)
                val in [".", "~"] && (hits += 1)
                print(io, val)
            end
            print(String(take!(io)))
            print(" h=$hits n = $(sum(calls)) max=$(maximum(calls))")
            println("")
        end
    end
end


