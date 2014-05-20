## Find real zeros of a polynomial
## This follows simplest algorithm outlined here: http://en.wikipedia.org/wiki/Vincent's_theorem
## Basically there are these steps:
## * we make the polynomial square free by finding psf = gcd(p, p') using gcd code from multroot
## * we consider the positive roots of the polynomial x->psf(x) and x->psf(-x)
## * for each an upper bound on the maximal real root is found using the Cauchy method.
##   other tighter bounds are possible, but this is easier to implement
## * we then do a bisection method on [0, ub]. The number of real roots in [a,b] is bounded above
##   by descarte's sign change rule on a modified polynomial. It this is 0, there is no root, 1 there is one
##   root we find using fzero, more than 1, then we bisect interval. 
##
## Drawbacks: the finding of a gcd to make p square free is problematic:
## * It can fail for large polynomials (wilkinson)
## * It can fail for nearby roots. (x-1+delta)*(x-1-delta) type things

## Some utilities

## hacky way to get coefficients from a power series object [a0 a1 ... an] -> a0 + a1x + a2x^2 + ... + an x^n
function coefficients(p::PowerSeries.AbstractSeries)
    n = length(names(p))
    map(i -> p.(symbol(string("c", i))), 0:(n-1))
end


## degree of powerseries
function degree(p::PowerSeries.AbstractSeries)
    as = coefficients(p)
    while (length(as) > 0) && (as[end] == 0.0)
        pop!(as)
    end
    length(as) - 1
end

## convert a polynomial to a power series
function poly_to_series(p::Poly) 
    Roots.degree(p) > 7 && PowerSeries.generate(Roots.degree(p))
    series(reverse(p.a)...)
end



## from base.@horner
function horner(x, ps) # ps = [a0, a1, a2, ..., an]
    ex = ps[end]
    for i in length(ps)-1:-1:1
        ex = ps[i] + x * ex
    end
    ex
end

horner(x, p::PowerSeries.AbstractSeries) = horner(x, coefficients(p))

## compose two series p(q(x))
function series_compose(p::PowerSeries.AbstractSeries, q::PowerSeries.AbstractSeries)
    horner(q, p)
end

## Polynomial pieces

## descarte's rule of signs
## counts a bound on number of real roots in (0,Inf) or (-Inf, 0)
## #of 0 <= roots = b - 2k, k=0, 1, ...
function descartes_bound(p::PowerSeries.AbstractSeries, dir=:positive)
    ps = reverse(coefficients(p))
    b = map(sign, ps)
    if dir != :positive
        b[1:2:length(b)] = -b[1:2:length(b)]
    end
    sum(map(*, b[1:end-1], b[2:end]) .< 0)        
end

## http://en.wikipedia.org/wiki/Vincent's_theorem
## gives bound on number of real roots in [a,b]
## We need to push through series here. Don't know how to do so via Polynomial type alone
function ab_test(p::PowerSeries.AbstractSeries, a, b)
    @assert a < b

    n = degree(p)
    x = series(tuple(0.0, 1.0, ntuple(n-1, x->0.0)...)...)

    ## the polynomial (1+x)^n * p((a+bx)/(1+x)) is considered for its posive roots
    (1+x)^n * series_compose(p, (a + b*x)/(1+x)) |> descartes_bound
end
    
      
## cauchy upperbound on positive real roots
## there are tighter bounds...
## cf. http://www.jucs.org/jucs_13_4/a_comparison_of_various/jucs_13_4_0455_0467_akritas.pdf
##     http://ssmr.ro/bulletin/pdf/53-3/akritas.pdf
function upperbound(p::PowerSeries.AbstractSeries)
    as = coefficients(p)             # a0, a1, ..., an
    if as[end] < 0
        as = -as
    end
    n = length(as) - 1
    ps = map(sign, as)         # signs
    ks = ps .< 0               # negative coefficients
    lamda = sum(ks)
    
    lamda > 0 || return(0) # nothing to be found if no negative coefficients

    mx = 0.0
    for k in 1:n
        if as[n+1-k] < 0
            mx = max(mx, (-lamda * as[n+1-k]/ as[n+1])^(1/k))
        end
    end
    mx + 1e-2 * rand()          # make slightly large
end

    

## use bisection method -- easiest to implement
## positive roots are all in (0, ub)
## p is square-free
function VAG(p::PowerSeries.AbstractSeries, a, b) 
    if a > b
        a,b = b, a
    end

    rts = Float64[]
    nr = ab_test(p, a, b)

    if nr == 1
        push!(rts, fzero(x -> horner(x, p), [a, b]))
    elseif nr > 1
        gamma = (a + b)/2
        tmp = horner(gamma, p) 
        if horner(gamma, p) == 0.0
            append!(rts, [gamma])
        end
        if a< gamma < b
            append!(rts, VAG(p, a, gamma))
            append!(rts, VAG(p, gamma, b))
        end
    end
    
    return(rts)
end
    
    

## find values on left and right of 0
function VAG(p::PowerSeries.AbstractSeries; tol::Real=eps())
    rts = Float64[]

    if abs(horner(0.0, p)) <= tol
        push!(rts, 0.0)
    end

    ## positive roots
    append!(rts, VAG(p, 2*eps(), upperbound(p)))
    ## negative roots. Need p->-p
    as = coefficients(p)
    as[2:2:length(as)] = -as[2:2:length(as)]
    q = series(as...)
    append!(rts, -VAG(q, 2*eps(), upperbound(q)))
end
    
## find real_roots of a polynomial
function real_roots(p::Poly)
    try
        ## need square free
        n, u, v, w = Roots.gcd_degree(p)
        VAG(poly_to_series(v))
    catch e
        ## guess that thing is square free!
        VAG(poly_to_series(p))
    end
end
    

## tests
if false
    x = poly([0.0])             # x-0 polynomial

    ## multiple roots
    real_roots(x^10*(x-1)^10 * (x-2)^10)

    ## large degree
    n = 6; real_roots(poly([1.0:n]))
    n = 7; real_roots(poly([1.0:n])) # fails -- gcd code is poor
    
    ## small separatoin
    delta = 1e-3
    real_roots(x*(x-delta)*(x+delta))

    delta = 1e-4                # misses 0
    real_roots(x*(x-delta)*(x+delta))

end
