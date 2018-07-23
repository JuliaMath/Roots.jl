# some simpler (and faster) implementations for root finding
#
# Not exported
#
# These avoid the setup costs of the `find_zero` method, so should be faster
# though they will take similar number of function calls.
#
# `Roots.bisection(f, a, b)`  (Bisection). 
# `Roots.a42(f, a, b)` (Roots.A42) Alefeld, Potra & Shi algorithm, like Brent's 
# `Roots.secant_method(f, xs)` (Order1) secant method
# `Roots.dfree(f, xs)`  (Order0) more robust secant method
#

"""
    bisection(f, a, b; [xatol, xrtol])

Performs bisection method to find a zero of a continuous
function. 

It is assumed that (a,b) is a bracket, that is, the function has
different signs at a and b. The interval (a,b) is converted to floating point
and shrunk when a or b is infinite. The function f may be infinite for
the typical case. If f is not continuous, the algorithm may find
jumping points over the x axis, not just zeros.

When tolerances are specified, this can use far fewer function calls
than the `Bisection()` method for `find_zero`.
    
If non-trivial tolerances are specified, the process will terminate
when the distance between the midpoint and the current bracket exceeds
the combined tolerance `max(xatol, max(|a|,|b|)*xrtol)`.
    
"""
function bisection(f, a, b; xatol=zero(float(a)), xrtol=zero(one(float(a))))
    u, v = promote(float(a), float(b))
    if u > v
        u,v=v,u
    end
    
    if isinf(u)
        u = nextfloat(u)
    end
    
    if isinf(v)
        v = prevfloat(v)
    end
    
    sign(f(u)) * sign(f(v)) < 0 || throw(ArgumentError("Interval [a,b] is not a bracket. f(a) and f(b) must have opposite signs"))

    zero_tolerance =  iszero(max(xatol, oneunit(a) * xrtol))
    _bisect(Val{zero_tolerance},  f, u, v, xatol, xrtol)
end

# if zero tolerance send to bisection64 or a42
function _bisect(::Type{Val{true}}, f, a::T, b::T, xatol, xrtol) where {T <: FloatNN}
    bisection64(f, a, b)
end
function _bisect(::Type{Val{true}}, f, a::T, b::T, xatol, xrtol) where {T}
    a42(f, a, b)
end
# if non zero tolerance, use faster middle.
function _bisect(::Type{Val{false}},  f, a::T, b::T, xatol, xrtol) where {T}
    
    fa, fb = sign(f(a)), sign(f(b))
    
    iszero(fa) && return a
    iszero(fb) && return b


    m::T = a + 0.5 * (b-a)
    tol = max(xatol, max(abs(a), abs(b)) * xrtol)
    
    while a < m < b

        (m-a) < tol && break 
        (b-m) < tol && break 
        
        
        fm = sign(f(m))

        iszero(fm) && break
        isnan(fm) && break
        
        if fa * fm < 0
            b,fb = m, fm
        else
            a,fa = m, fm
        end

        m = a + 0.5*(b-a)
        
    end

    m
end

"""
    secant_method(f, xs; [atol=0.0, rtol=8eps(), maxevals=1000])

Perform secant method to solve f(x) = 0.

The secant method is an iterative method with update step
given by b - fb/m where m is the slope of the secant line between
(a,fa) and (b,fb).

The inital values can be specified as a pair of 2, as in `(a,b)` or
`[a,b]`, or as a single value, in which case a value of `b` is chosen.

The algorithm returns m when `abs(fm) <= max(atol, abs(m) * rtol)`.
If this doesn't occur before `maxevals` steps or the algorithm encounters
an issue, a value of `NaN` is returned

The `Order1` method for `find_zero` also implements the secant
method. This one will be faster, as there are fewer setup costs.

Examples
```julia
Roots.secant_method(sin, (3,4))
Roots.secant_method(x -> x^5 -x - 1, 1.1)
```
"""
function secant_method(f, xs; atol=zero(float(real(first(xs)))), rtol=8eps(one(float(real(first(xs))))), maxevals=100)

    if length(xs) == 1 # secant needs x0, x1; only x0 given
        a = float(xs[1])
        fa = f(a)

        h = eps(one(real(a)))^(1/3)
        da = h*oneunit(a) + abs(a)*h^2 # adjust for if eps(a) > h
        b = a + da
        
    else
        a, b = promote(float(xs[1]), float(xs[2]))
    end
    secant(f, a, b, atol, rtol, maxevals)
end


function secant(f, a::T, b::T, atol=zero(T), rtol=8eps(T), maxevals=100) where {T}
    nan = (0a)/(0a)
    cnt = 0

    fa, fb = f(a), f(b)
    fb == fa && return nan

    uatol = atol / oneunit(atol) * oneunit(real(a))
    adjustunit = oneunit(real(fb))/oneunit(real(b))

    while cnt < maxevals
        m = b - (b-a)*fb/(fb - fa)
        fm = f(m)

        iszero(fm) && return m
        isnan(fm) && return nan
        abs(fm) <= adjustunit * max(uatol, abs(m) * rtol) && return m
        if fm == fb
            sign(fm) * sign(f(nextfloat(m))) <= 0 && return m
            sign(fm) * sign(f(prevfloat(m))) <= 0 && return m            
            return nan
        end
        
        a,b,fa,fb = b,m,fb,fm

        cnt += 1
    end

    return nan
end 



## This is basically Order0(), but with different, default, tolerances employed
## It takes more function calls, but works harder to find exact zeros
## where exact means either iszero(fx), adjacent floats have sign change, or
## abs(fxn) <= 8 eps(xn)
"""
    dfree(f, xs)

A more robust secant method implementation

Solve for `f(x) = 0` using an alogorithm from *Personal Calculator Has Key
to Solve Any Equation f(x) = 0*, the SOLVE button from the
[HP-34C](http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf).

This is also implemented as the `Order0` method for `find_zero`.

The inital values can be specified as a pair of two values, as in
`(a,b)` or `[a,b]`, or as a single value, in which case a value of `b`
is computed, possibly from `fb`.  The basic idea is to follow the
secant method to convergence unless:

* a bracket is found, in which case bisection is used;

* the secant method is not converging, in which case a few steps of a
  quadratic method are used to see if that improves matters.

Convergence occurs when `f(m) == 0`, there is a sign change between
`m` and an adjacent floating point value, or `f(m) <= 2^3*eps(m)`.
      
A value of `NaN` is returned if the algorithm takes too many steps
before identifying a zero.

# Examples

```julia
Roots.dfree(x -> x^5 - x - 1, 1.0)
```

"""
function dfree(f, xs)

    if length(xs) == 1
        a = float(xs[1])
        fa = f(a)

        h = eps(one(a))^(1/3)
        da = h*oneunit(a) + abs(a)*h^2 # adjust for if eps(a) > h
        b = float(a + da)
        fb = f(b)
    else
        a, b = promote(float(xs[1]), float(xs[2]))
        fa, fb = f(a), f(b)
    end

    
    nan = (0*a)/(0*a) # try to preserve type
    cnt, MAXCNT = 0, 5 * ceil(Int, -log(eps(one(a))))  # must be higher for BigFloat
    MAXQUAD = 3
    
    if abs(fa) > abs(fb)
        a,fa,b,fb=b,fb,a,fa
    end
    
    # we keep a, b, fa, fb, gamma, fgamma
    quad_ctr = 0
    while !iszero(fb)
        cnt += 1

        if sign(fa) * sign(fb) < 0
            return bisection(f, a, b)
        end

        # take a secant step
        gamma =  float(b - (b-a) * fb / (fb - fa))
        # modify if gamma is too small or too big
        if iszero(abs(gamma-b))
            gamma = b + 1/1000 * abs(b-a)  # too small
        elseif abs(gamma-b)  >= 100 * abs(b-a)
            gamma = b + sign(gamma-b) * 100 * abs(b-a)  ## too big
        end
        fgamma = f(gamma)
        
        # change sign
        if sign(fgamma) * sign(fb) < 0
            return bisection(f, gamma, b)
        end 

        # decreasing
        if abs(fgamma) < abs(fb)
            a,fa, b,fb = b, fb, gamma, fgamma
            quad_ctr = 0
            cnt < MAXCNT && continue
        end

        gamma = float(quad_vertex(a,fa,b,fb,gamma,fgamma))
        fgamma = f(gamma)
        # decreasing now?
        if abs(fgamma) < abs(fb)
            a,fa, b,fb = b, fb, gamma, fgamma
            quad_ctr = 0
            cnt < MAXCNT && continue
        end

        
        quad_ctr += 1
        if (quad_ctr > MAXQUAD) || (cnt > MAXCNT) || iszero(gamma - b)  || isnan(gamma)
            bprev, bnext = prevfloat(b), nextfloat(b)
            fbprev, fbnext = f(bprev), f(bnext)
            sign(fb) * sign(fbprev) < 0 && return b
            sign(fb) * sign(fbnext) < 0 && return b
            for (u,fu) in ((b,fb), (bprev, fbprev), (bnext, fbnext))
                abs(fu)/oneunit(fu) <= 2^3*eps(u/oneunit(u)) && return u
            end
            return nan # Failed. 
        end
        
        if abs(fgamma) < abs(fb)
            b,fb, a,fa = gamma, fgamma, b, fb
        else
            a, fa = gamma, fgamma
        end
        
    end
    b

end

