###


const bracketing_error = """The interval [a,b] is not a bracketing interval.
You need f(a) and f(b) to have different signs (f(a) * f(b) < 0).
Consider a different bracket or try fzero(f, c) with an initial guess c.

"""

## Methods for root finding which use a bracket

## Bisection for FLoat64 values.
##
## From Jason Merrill https://gist.github.com/jwmerrill/9012954
## cf. http://squishythinking.com/2014/02/22/bisecting-floats/
# Alternative "mean" definition that operates on the binary representation
# of a float. Using this definition, bisection will never take more than
# 64 steps (over Float64)
const FloatNN = Union{Float64, Float32, Float16}
const _float_int_pairs = Dict(Float64 => UInt64, Float32 => UInt32, Float16 => UInt16)

function _middle(x::T, y::T) where {T <: FloatNN}
    # Use the usual float rules for combining non-finite numbers
    if !isfinite(x) || !isfinite(y)
        return x + y
    end
    # Always return 0.0 when inputs have opposite sign
    if sign(x) != sign(y) && !iszero(x) && !iszero(y)
        return zero(T)
    end
    
    negate = sign(x) < 0 || sign(y) < 0 

    # do division over unsigned integers with bit shift
    xint = reinterpret(_float_int_pairs[T], abs(x))
    yint = reinterpret(_float_int_pairs[T], abs(y))
    mid = (xint + yint) >> 1

    # reinterpret in original floating point
    unsigned = reinterpret(T, mid)

    negate ? -unsigned : unsigned
end

## fall back or non Floats
function _middle(x::Number, y::Number)
    x + (y-x)/2
end


"""

    Roots.bisection64(f, a, b)

(unexported)

* `f`: a callable object, like a function

* `a`, `b`: Real values specifying a *bracketing* interval (one with
`f(a) * f(b) < 0`). These will be converted to `Float64` values.

Runs the bisection method using midpoints determined by a trick
leveraging 64-bit floating point numbers. After ensuring the
intermediate bracketing interval does not straddle 0, the "midpoint"
is half way between the two values onces converted to unsigned 64-bit
integers. This means no more than 64 steps will be taken, fewer if `a`
and `b` already share some bits.

The process is guaranteed to return a value `c` with `f(c)` one of
`0`, `Inf`, or `NaN`; *or* one of `f(prevfloat(c))*f(c) < 0` or
`f(c)*f(nextfloat(c)) > 0` holding. 

This function is a bit faster than the slightly more general 
`find_zero(f, [a,b], Bisection())` call.

Due to Jason Merrill.

"""
function bisection64(@nospecialize(f), a0::FloatNN, b0::FloatNN)

    a,b = promote(a0, b0)


    if a > b
        b,a = a, b
    end
    
    
    m = _middle(a,b)

    fa, fb = sign(f(a)), sign(f(b))

    
    fa * fb > 0 && throw(ArgumentError(bracketing_error)) 
    (iszero(fa) || isnan(fa) || isinf(fa)) && return a
    (iszero(fb) || isnan(fb) || isinf(fb)) && return b
    
    while a < m < b

        fm = sign(f(m))

        if iszero(fm) || isnan(fm) || isinf(fm)
            return m
        elseif fa * fm < 0
            b,fb=m,fm
        else
            a,fa=m,fm
        end
        m = _middle(a,b)
    end
    return m
end


####
## find_zero interface.
"""

    Bisection()

If possible, will use the bisection method over `Float64` values. The
bisection method starts with a bracketing interval `[a,b]` and splits
it into two intervals `[a,c]` and `[c,b]`, If `c` is not a zero, then
one of these two will be a bracketing interval and the process
continues. The computation of `c` is done by `_middle`, which
reinterprets floating point values as unsigned integers and splits
there. This method avoids floating point issues and guarantees a
"best" solution (one where a zero is found or the bracketing interval
is of the type `[a, nextfloat(a)]`).

When this is not possible, will default to `A42` method.
    
"""    
mutable struct Bisection <: AbstractBisection end
mutable struct Bisection64 <: AbstractBisection end

"""
    Roots.A42()

Bracketing method which finds the root of a continuous function within
a provided interval [a, b], without requiring derivatives. It is based
on algorithm 4.2 described in: 1. G. E. Alefeld, F. A. Potra, and
Y. Shi, "Algorithm 748: enclosing zeros of continuous functions," ACM
Trans. Math. Softw. 21, 327–344 (1995).
"""
mutable struct A42 <: AbstractBisection end


function init_options(::M,
                      state;
                      xatol=missing,
                      xrtol=missing,
                      atol=missing,
                      rtol=missing,
                      maxevals::Int=typemax(Int),
                      maxfnevals::Int=typemax(Int),
                      verbose::Bool=false,
                      kwargs...) where {M <: Union{Bisection64, A42}}

    ## Where we set defaults
    x1 = real(oneunit(state.xn1))
    fx1 = real(oneunit(float(state.fxn1)))

    ## map old tol names to new
    ## deprecate in future
    xatol, xrtol, atol, rtol = _map_tolerance_arguments(Dict(kwargs), xatol, xrtol, atol, rtol)
    
    # all are 0 by default
    options = UnivariateZeroOptions(ismissing(xatol) ? zero(x1) : xatol,       # unit of x
                                    ismissing(xrtol) ?  zero(x1/oneunit(x1)) : xrtol,               # unitless
                                    ismissing(atol)  ? zero(fx1) : atol,  # units of f(x)
                                    ismissing(rtol)  ?  zero(fx1/oneunit(fx1)) : rtol,            # unitless
                                    maxevals, maxfnevals,
    verbose)    

    options
end

## we dispatch to either floating-point-bisection or A42 here.
function find_zero(fs, x0, method::AbstractBisection; kwargs...)
    
    x = float.(x0)
    F = callable_function(fs)
    state = init_state(method, F, x)
    options = init_options(method, state; kwargs...)

    # we try a faster alternative for floating point values based on verboseness
    isa(method, A42) && return find_zero(method, F, options, state)
    isa(method, FalsePosition) && return find_zero(method, F, options, state)
    
    T = eltype(state.xn1)
    if T <: FloatNN
        if options.verbose
            find_zero(Bisection64(), F, options, state)
        else
            x0, x1 = state.xn0, state.xn1
            state.xn1 = bisection64(F, x0, x1)
            state.message = "Used bisection to find the zero, steps not counted."
            state.stopped = state.x_converged = true
            return state.xn1
        end
    else
        find_zero(A42(), F, options, state)
    end
    
end


function find_zero(method::A42, F, options::UnivariateZeroOptions, state::UnivariateZeroState{T,S}) where {T<:Number, S<:Number}
     x0, x1 = state.xn0, state.xn1
    state.xn1 = a42(F, x0, x1; xtol=options.xabstol, maxeval=options.maxevals,
                        verbose=options.verbose)
        state.message = "Used Alefeld-Potra-Shi method, `Roots.a42`, to find the zero. Iterations and function evaluations are not counted properly."
        state.stopped = state.x_converged = true
        
        options.verbose && show_trace(state, [state.xn1], [state.fxn1], method)
    return state.xn1
end


## in Order0, we run bisection if a bracketing interval is found
## this is meant to be as speedy as possible
function _run_bisection(fs, options, state::UnivariateZeroState{T,S}) where {T<:FloatNN, S<:Number}
    xn0, xn1 = state.xn0, state.xn1
    state.xn1 = bisection64(fs, xn0, xn1)
    state.x_converged = true
    state.message = "Used Bisection() for last step, steps not counted"
end

function _run_bisection(fs, options, state::UnivariateZeroState{T,S}) where {T<:Number, S<:Number}
    state.xn1 = find_zero(A42(), fs, options, state)
    state.x_converged = true
    state.message = "Used A42() for last step, steps not counted"
end


# ## helper function
function adjust_bracket(x0)
    u, v = float.(promote(x0...))
    if u > v
        u, v = v, u
    end


    if isinf(u)
        u = nextfloat(u)
    end
    if isinf(v)
        v = prevfloat(v)
    end
    u, v
end



function init_state(method::AbstractBisection, fs, x)
    length(x) > 1 || throw(ArgumentError(bracketing_error))
    
    x0, x2 = adjust_bracket(x)
    y0, y2 = promote(fs(x0), fs(x2))

    sign(y0) * sign(y2) > 0 && throw(ArgumentError(bracketing_error))

    state = UnivariateZeroState(x0, x2,
                                y0, y2,
                                0, 2,
                                false, false, false, false,
                                "")
    state
end


## This uses _middle bisection Find zero using modified bisection
## method for FloatXX arguments.  This is guaranteed to take no more
## steps the bits of the type. The `a42` alternative usually has fewer
## iterations, but this seems to find the value with fewer function
## evaluations.
##
## This terminates when there is no more subdivision or function is zero

function update_state(method::Union{Bisection,Bisection64}, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T<:Number,S<:Number} 
    x0, x2 = o.xn0, o.xn1
    y0, y2 = o.fxn0, o.fxn1

    x1 = _middle(x0, x2)

    y1 = fs(x1)
    incfn(o)

    if sign(y0) * sign(y1) > 0
        x0, x2 = x1, x2
        y0, y2 = y1, y2
    else
        x0, x2 = x0, x1
        y0, y2 = y0, y1
    end

    o.xn0, o.xn1 = x0, x2
    o.fxn0, o.fxn1 = y0, y2
    incsteps(o)
    nothing
end

## convergence is much different here
## the method converges,
## as we bound between x0, nextfloat(x0) is not measured by eps(), but eps(x0)
function assess_convergence(method::Union{Bisection64,Bisection}, state, options)
    x0, x2 = state.xn0, state.xn1
    if x0 > x2
        x0, x2 = x2, x0
    end

    x1 = _middle(x0, x2)
    if iszero(state.fxn1)
        state.message = ""
        state.stopped = state.f_converged = true
        return true
    end

    x0 == x1 || x1 == x2 && return true
    d1 = isapprox(x0, x1, atol=options.xabstol, rtol=options.xreltol)
    d2 = isapprox(x1, x2, atol=options.xabstol, rtol=options.xreltol)
    !d1 && !d2 && return false
    
#    x0 < x1 && x1 < x2 && return false
     

    state.message = ""
    state.stopped = state.x_converged = true
    true
end


##################################################

"""

    Roots.a42(f, a, b; kwargs...)

(not exported)

Finds the root of a continuous function within a provided
interval [a, b], without requiring derivatives. It is based on algorithm 4.2
described in: 1. G. E. Alefeld, F. A. Potra, and Y. Shi, "Algorithm 748:
enclosing zeros of continuous functions," ACM Trans. Math. Softw. 21,
327–344 (1995).


input:
    `f`: function to find the root of
    `a`, `b`: the initial bracket, with: a < b, f(a)*f(b) < 0
    `xtol`: acceptable error (it's safe to set zero for machine precision)
    `maxeval`:  maximum number of iterations

output:
    an estimate of the zero of f

By John Travers

"""
function a42(f, a, b;
      xtol=zero(float(a)),
      maxeval::Int=15,
      verbose::Bool=false)

    if a > b
        a,b = b,a
    end
    fa, fb = f(a), f(b)

    if a >= b || sign(fa)*sign(fb) >= 0
        error("on input a < b and f(a)f(b) < 0 must both hold")
    end
    if xtol/oneunit(xtol) < 0.0
        error("tolerance must be >= 0.0")
    end
    
    c, fc = secant(f, a, fa, b, fb)
    a42a(f, float(a), fa, float(b), fb, float(c), fc,
         xtol=xtol, maxeval=maxeval, verbose=verbose)
end

"""

Split Alefeld, F. A. Potra, and Y. Shi algorithm 4.2 into a function
where `c` is passed in.

Solve f(x) = 0 over bracketing interval [a,b] starting at c, with a < c < b

"""
a42a(f, a, b, c=(a+b)/2; args...) = a42a(f, a, f(a), b, f(b), c, f(c); args...)

function a42a(f, a, fa, b, fb, c, fc;
       xtol=zero(float(a)),
       maxeval::Int=15,
       verbose::Bool=false)

    try
        # re-bracket and check termination
        a, fa, b, fb, d, fd = bracket(f, a, fa, b, fb, c, fc, xtol)
        ee, fee = d, fd
        for n = 2:maxeval
            # use either a cubic (if possible) or quadratic interpolation
            if n > 2 && distinct(a, fa, b, fb, d, fd, ee, fee)
                c, fc = ipzero(f, a, fa, b, fb, d, fd, ee, fee)
            else
                c, fc = newton_quadratic(f, a, fa, b, fb, d, fd, 2)
            end
            # re-bracket and check termination
            ab, fab, bb, fbb, db, fdb = bracket(f, a, fa, b, fb, c, fc, xtol)
            eb, feb = d, fd
            # use another cubic (if possible) or quadratic interpolation
            if distinct(ab, fab, bb, fbb, db, fdb, eb, feb)
                cb, fcb = ipzero(f, ab, fab, bb, fbb, db, fdb, eb, feb)
            else
                cb, fcb = newton_quadratic(f, ab, fab, bb, fbb, db, fdb, 3)
            end
            # re-bracket and check termination
            ab, fab, bb, fbb, db, fdb = bracket(f, ab, fab, bb, fbb, cb, fcb, xtol)
            # double length secant step; if we fail, use bisection
            if abs(fab) < abs(fbb)
                u, fu = ab, fab
            else
                u, fu = bb, fbb
            end
            # u = abs(fab) < abs(fbb) ? ab : bb
#            cb = u - 2*fu/(fbb - fab)*(bb - ab)
            del =  2*fu/(fbb - fab)*(bb - ab)
            if !isnan(del) # add check on NaN
                cb = u - del
            end
            fcb = f(cb)
            if abs(cb - u) > (bb - ab)/2
                ch, fch = ab+(bb-ab)/2, f(ab+(bb-ab)/2)
            else
                ch, fch = cb, fcb
            end
            # ch = abs(cb - u) > (bb - ab)/2 ? ab + (bb - ab)/2 : cb
            # re-bracket and check termination
            ah, fah, bh, fbh, dh, fdh = bracket(f, ab, fab, bb, fbb, ch, fch, xtol)
            # if not converging fast enough bracket again on a bisection
            if bh - ah < 0.5*(b - a)
                a, fa = ah, fah
                b, fb = bh, fbh
                d, fd = dh, fdh
                ee, fee = db, fdb
            else
                ee, fee = dh, fdh
                a, fa, b, fb, d, fd = bracket(f, ah, fah, bh, fbh, ah + (bh - ah)/2, f(ah+(bh-ah)/2), xtol)
            end

            verbose && println(@sprintf("a=%18.15f, n=%s", float(a), n))

            if nextfloat(float(ch)) * prevfloat(float(ch)) <= 0 * oneunit(ch)^2
                throw(StateConverged(ch))
            end
            if nextfloat(float(a)) >= b
                throw(StateConverged(a))
            end
        end
        throw(ConvergenceFailed("More than $maxeval iterations before convergence"))
    catch ex
        if isa(ex, StateConverged)
            return ex.x0
        else
            rethrow(ex)
        end
    end
    throw(ConvergenceFailed("More than $maxeval iterations before convergence"))
end



# calculate a scaled tolerance
# based on algorithm on page 340 of [1]
function tole(a::S, b::R, fa, fb, tol) where {S,R}
    u = abs(fa) < abs(fb) ? abs(a) : abs(b)
    T = promote_type(S,R)
    epsilon = eps(T)
    2u*(epsilon/oneunit(epsilon)) + tol
end


# bracket the root
# inputs:
#     - f: the function
#     - a, b: the current bracket with a < b, f(a)f(b) < 0
#     - c within (a,b): current best guess of the root
#     - tol: desired accuracy
#
# if root is not yet found, return
#     ab, bb, d
# with:
#     - [ab, bb] a new interval within [a, b] with f(ab)f(bb) < 0
#     - d a point not inside [ab, bb]; if d < ab, then f(ab)f(d) > 0,
#       and f(d)f(bb) > 0 otherwise
#
# if the root is found, throws a StateConverged instance with x0 set to the
# root.
#
# based on algorithm on page 341 of [1]
function bracket(f, a, fa, b, fb, c, fc, tol)
    
    if !(a <= c <= b)
        error("c must be in (a,b)")
    end
    delta = 0.7*tole(a, b, fa, fb, tol)
    if b - a <= 4delta
        c = (a + b)/2
        fc = f(c)
    elseif c <= a + 2delta
        c = a + 2delta
        fc = f(c)
    elseif c >= b - 2delta
        c = b - 2delta
        fc = f(c)
    end
    if fc == 0
        throw(Roots.StateConverged(c))
    elseif sign(fa)*sign(fc) < 0 
        aa, faa = a, fa
        bb, fbb = c, fc
        db, fdb = b, fb
    else
        aa, faa = c, fc
        bb, fbb = b, fb
        db, fdb = a, fa
    end
    if bb - aa < 2*tole(aa, bb, faa, fbb, tol)
        x0 = abs(faa) < abs(fbb) ? aa : bb
        throw(Roots.StateConverged(x0))
    end
    aa, faa, bb, fbb, db, fdb
end


# take a secant step, if the resulting guess is very close to a or b, then
# use bisection instead
function secant(f, a::T, fa, b, fb) where {T}

    c = a - fa/(fb - fa)*(b - a)
    _tol = 5*eps(T)
    tol=_tol/oneunit(_tol)

    if isnan(c) || c <= a + abs(a)*tol || c >= b - abs(b)*tol
        return a + (b - a)/2, f(a+(b-a)/2)
    end
    return c, f(c)
end


# approximate zero of f using quadratic interpolation
# if the new guess is outside [a, b] we use a secant step instead
# based on algorithm on page 330 of [1]
function newton_quadratic(f, a, fa, b, fb, d, fd, k::Int)
    B = (fb - fa)/(b - a)
    A = ((fd - fb)/(d - b) - B)/(d - a)
    if A == 0
        return secant(f, a, fa, b, fb)
    end
    
    r = A*fa/oneunit(A*fa) > 0 ? a : b
    for i = 1:k
        r -= (fa + (B + A*(r - b))*(r - a))/(B + A*(2*r - a - b))
    end
    if isnan(r) || (r <= a || r >= b)
        r, fr = secant(f, a, fa, b, fb)
    else
        fr = f(r)
    end
    return r, fr
end


# approximate zero of f using inverse cubic interpolation
# if the new guess is outside [a, b] we use a quadratic step instead
# based on algorithm on page 333 of [1]
function ipzero(f, a, fa, b, fb, c, fc, d, fd)

    Q11 = (c - d)*fc/(fd - fc)
    Q21 = (b - c)*fb/(fc - fb)
    Q31 = (a - b)*fa/(fb - fa)
    D21 = (b - c)*fc/(fc - fb)
    D31 = (a - b)*fb/(fb - fa)
    Q22 = (D21 - Q11)*fb/(fd - fb)
    Q32 = (D31 - Q21)*fa/(fc - fa)
    D32 = (D31 - Q21)*fc/(fc - fa)
    Q33 = (D32 - Q22)*fa/(fd - fa)
    c = a + (Q31 + Q32 + Q33)
    if (c <= a) || (c >= b)
        return newton_quadratic(f, a, fa, b, fb, d, fd, 3)
    end
    return c, f(c)
end


# floating point comparison function
function almost_equal(x::T, y::T) where {T}
    min_diff = oneunit(x) * realmin(float(x/oneunit(x)))*32
    abs(x - y) < min_diff
end


# check that all interpolation values are distinct
function distinct(a, f1, b, f2, d, f3, e, f4)
    !(almost_equal(f1, f2) || almost_equal(f1, f3) || almost_equal(f1, f4) ||
      almost_equal(f2, f3) || almost_equal(f2, f4) || almost_equal(f3, f4))
end

## ----------------------------

"""

    FalsePosition()

Use the [false
position](https://en.wikipedia.org/wiki/False_position_method) method
to find a zero for the function `f` within the bracketing interval
`[a,b]`.

The false position method is a modified bisection method, where the
midpoint between `[a_k, b_k]` is chosen to be the intersection point
of the secant line with the x axis, and not the average between the
two values.

To speed up convergence for concave functions, this algorithm
implements the 12 reduction factors of Galdino (*A family of regula
falsi root-finding methods*). These are specified by number, as in
`FalsePosition(2)` or by one of three names `FalsePosition(:pegasus)`,
`FalsePosition(:illinois)`, or `FalsePosition(:anderson_bjork)` (the
default). The default choice has generally better performance than the
others, though there are exceptions.

For some problems, the number of function calls can be greater than
for the `bisection64` method, but generally this algorithm will make
fewer function calls.

Examples
```
find_zero(x -> x^5 - x - 1, [-2, 2], FalsePosition())
```
"""    
struct FalsePosition <: AbstractBisection
    reduction_factor::Union{Int, Symbol}
    FalsePosition(x=:anderson_bjork) = new(x)
end

function update_state(method::FalsePosition, fs, o, options)

    fs
    a, b =  o.xn0, o.xn1

    fa, fb = o.fxn0, o.fxn1

    lambda = fb / (fb - fa)
    tau = 1e-10                   # some engineering to avoid short moves
    if !(tau < abs(lambda) < 1-tau)
        lambda = 1/2
    end
    x = b - lambda * (b-a)        
    fx = fs(x)
    incfn(o)
    incsteps(o)

    if iszero(fx)
        o.xn1 = x
        o.fxn1 = fx
        return
    end

    if sign(fx)*sign(fb) < 0
        a, fa = b, fb
    else
        fa = galdino[method.reduction_factor](fa, fb, fx)
    end
    b, fb = x, fx

    
    o.xn0, o.xn1 = a, b 
    o.fxn0, o.fxn1 = fa, fb
    
    nothing
end

# the 12 reduction factors offered by Galadino
galdino = Dict{Union{Int,Symbol},Function}(:1 => (fa, fb, fx) -> fa*fb/(fb+fx),
                                           :2 => (fa, fb, fx) -> (fa - fb)/2,
                                           :3 => (fa, fb, fx) -> (fa - fx)/(2 + fx/fb),
                                           :4 => (fa, fb, fx) -> (fa - fx)/(1 + fx/fb)^2,
                                           :5 => (fa, fb, fx) -> (fa -fx)/(1.5 + fx/fb)^2,
                                           :6 => (fa, fb, fx) -> (fa - fx)/(2 + fx/fb)^2,
                                           :7 => (fa, fb, fx) -> (fa + fx)/(2 + fx/fb)^2,
                                           :8 => (fa, fb, fx) -> fa/2,
                                           :9 => (fa, fb, fx) -> fa/(1 + fx/fb)^2,
                                           :10 => (fa, fb, fx) -> (fa-fx)/4,
                                           :11 => (fa, fb, fx) -> fx*fa/(fb+fx),
                                           :12 => (fa, fb, fx) -> (fa * (1-fx/fb > 0 ? 1-fx/fb : 1/2))  
)
# give common names
for (nm, i) in [(:pegasus, 1), (:illinois, 8), (:anderson_bjork, 12)]
    galdino[nm] = galdino[i]
end

## --------------------------------------

## find zeros code
## Interval is either a bracket or not
struct Interval
l
m
r
fl
fm
fr
end

midpoint(l, r) = l + 0.484375 * (r-l) # avoid 1/2 to avoid nice fractional power
function Interval(f, l, r)
    m = midpoint(l, r) 
    Interval(l, m, r, f(l), f(m), f(r))
end

function split_interval(f, a)
    lm = midpoint(a.l, a.m)
    rm = midpoint(a.m, a.r)
    flm, frm = f(lm), f(rm)
    Interval(a.l, lm, a.m, a.fl, flm, a.fm), Interval(a.m, rm, a.r, a.fm, frm, a.fr)
end

isbracket(a::Interval) = sign(a.fl) * sign(a.fr) < 0
deltax(a::Interval) = a.r - a.l
deltaf(a::Interval) = a.fr - a.fl

function first_derivative(a::Interval)
    deltaf(a) / deltax(a)
end

# from adaptive_grid
function second_derivative(a::Interval)
    x1,x2,x3 = a.l, a.m, a.r
    y1,y2,y3 = a.fl, a.fm, a.fr
    fpp = 2 * (y1/(x2-x1)/(x3-x1) - y2/(x3-x2)/(x2-x1) + y3/(x3-x2)/(x3-x1))
    fpp
end


# Halving threshold heuristic
# return true to keep, false to drop
function HT(a::Interval, C, maxmultiplicity=1)
    # why should we keep non-bracketing interval?
    # - quadratic expansion (with boost) is a bracket
    # - f is reasonably close
    dx = deltax(a) |> strip_unit
    fp = first_derivative(a) |> strip_unit
    fpp = second_derivative(a) |> strip_unit
    fmin = min(abs(a.fl), abs(a.fr)) |> strip_unit
    
    qstep = C * (abs(fp) * dx + 1/2 * abs(fpp) * dx^2)
    if isfinite(fpp) && fmin <= qstep 
        return true
    end
    
    # condition from paper
    # dx > 1/C * min(|fl|, |fr|) / dx
    lambda = min(1.0, abs(fpp)) 
    if fmin < lambda  * C * dx * dx^maxmultiplicity
        return true
    end

    return false
end

# check if a non-bracketing interval has any promise. Keep
# if it does. 
function has_promise(a::Interval, C, maxmultiplicity)
    HT(a,C, maxmultiplicity)
end

# have more data here, can better estimate?
function has_promise(a::Interval, la::Interval, ra::Interval, C, maxmultiplicity)
    HT(a, C, maxmultiplicity)
end
    
# bracket threshod
function BT(a::Interval, xatol, xrtol)
    # paper uses min(rtol * deltax(a), xatol), we use |a.m|
    # Here we want to ensure separation -- can bracket
    # have two or more distinct zeros

    max(xatol, abs(a.m) * xrtol)
end



## is fm ~ 0?
## A tricky question. In general, fm ~ 0 is handled by |fm| < atol, as rtol won't
## come into play due to |fm| factor.
## 
## However, if x is a root and x(1+d) the fp approx, we have the approximation
## f(x+xd) ~ f(x) + f'(x)*x*d, so if f(x) is zero
## and f(x(1+d)) the function eval on floating point approx, we have a tolerance
## that should be like C |x| * d. 
##
## If so, we then will use `find_zero` to search for a root in the interval
## As such, this uses a relaxed tolerance composed from rtol and atol
function isapproxzero(alpha::Interval, atol, rtol)

    # this seems wrong, but results are good on most test cases
#    return isapprox(alpha.fm, zero(alpha.fm), atol=atol, rtol=rtol)

    # this adds in factor of |x|, but keeps |fx| factor too
#    lambda = max(1.0, strip_unit(abs(alpha.m)))
#    return isapprox(alpha.fm, zero(alpha.fm), atol=atol, rtol=lambda * rtol)


    # this seems right but... not include |fx| with rtol
    axm, afm = abs(alpha.m), abs(alpha.fm)
    lambda = max(1.0, strip_unit(axm)) * oneunit(afm)
    ftol = max(atol, lambda * rtol)
    return afm < ftol

    # # more relaxed, lets find_zero polish off answer
    # lambda = max(1.0, strip_unit(axm)) * strip_unit(afm)
    # p = 1/2
    # ftol = lambda * rtol
    # ftol = (ftol < 1 ? (ftol)^(p) : ftol^(1/p))
    # tol = max(strip_unit(atol), ftol)
    # return strip_unit(afm) <= tol
end


"""
    
    find_zeros(f, a, b; C, maxmultiplicity, no_pts, xatl, xrtol, atol, rtol, maxevals, maxfnevals)

Search for zeros of a univariate function over the interval (a,b).

The interval (a,b) is broken into `no_pts` subintervals. Each interval
is split in two. Those which are bracketing intervals are kept, those
which show "promise" to contain a zero are kept, others are dropped
from consideration. The process repeats.

Once bracketing intervals are smaller than the proposed zero
separation, bisection is used to find the zero, guaranteed by the
bracket (assuming a continuous function, otherwise only a zero
crossing is guaranteed). Once "promising" intervals have a function
value at its midpoint close enough to a zero, `find_zero` will be used
to search for a zero in the interval.

The basic algorithm follows one proposed by
[Razbani](https://arxiv.org/pdf/1501.05298.pdf).

    
The adaptive selection of promising intervals is dependent on some parameters

* `xatol`, `xrtol`: These are used to set the size of the gap between
  successive zeros. (These are *not* used as a criteria for
  convergence, as when used with `find_zero`.) Bisection is used to
  find zeros of bracketing intervals once an interval is smaller than
  this gap. The number of function evaluations is dependent on this
  tolerance.


* `atol`, `rtol`: These are used to find a tolerance on the smallness
  of `f` at the midpoint of an interval. Once the value is smaller
  than the tolerance, `find_zero` is used to
  search for a zero within the interval.


* `C`, `maxmultiplicity`: An interval shows "promise" if an optimistic
  (quadratic) Taylor approximation indicates a possible zero *or* if
  `f` at the midpoint is smaller than ` lambda * C * dx *
  dx^maxmultiplicity` where `lambda` depends on an estimate for the
  second derivative. The heuristic is that the `f` value is reasonably
  small, so may lead to a zero, even if the interval is not a
  bracketing interval.      

The parameters `maxevals` and `maxfnevals` can be set to limit the
number of steps or function evaluations the algorithm will take.

Setting `verbose=true` will display a print out of the algorithm.


### Notes:

This zero finding algorithm employs heuristics that may not be
effective in certain cases. Function patterns that can cause these
heuristics to fail include:

* Non-simple zeros. A simple zero is one where f(x)/(x-c) does not
  have a zero at c. Non-simple zeros, such as `f(x) = x^2` may not
  have a bracketing and may be very close to 0 without having a roots,
  e.g. `f(x) = x^2 + 1e-10`. Increasing the `maxmultiplicity`
  parameter may help identify these zeros, though at the expense of
  more function evaluations.

* Close by zeros. The minimum gap between successive zeros, m, is an
  important criteria for a successful search. If this gap is too
  small, then a zero may be missed. In general, keeping (b-a)/m small
  helps the algorithm

This algorithm can require many function calls. Adjusting `C` to be
smaller can reduce this, but at the expense of more aggressive
definition of "promising." Similarly, reducing the ruler for the
smallest gap between successive zeros by adjusting `xatol` and `xrtol`
can reduce the number of function calls at the expense of not being
able to separate nearby zeros.

Finally, for some examples, small adjustments in the tolerances can lead quite
different results.
    
### Examples:

```
f(x) = exp(x) - x^4
find_zeros(f, -5, 20) # three zeros 8.61317, 1.42961, -0.815553
f(x) = cos(x) - x/10
find_zeros(f, 0, 10)  # three zeros 7.06889, 5.26712, 1.42755
```

    
"""    
function find_zeros(f, a, b;
                    C= 10 * one(float(a)),
                    maxmultiplicity = 1,
                    no_pts::Int=101,
                    xatol = 1e-3 * oneunit(float(a)),
                    xrtol = 1e-3,
                    atol = 8*eps(eltype(float(a))) * oneunit(f(float(a))),
                    rtol = eps(eltype(float(a))),
                    maxevals::Int=typemax(Int),
                    maxfnevals::Int=1_000_000,  # high, but not impossible to reach
                    verbose::Bool=false
    )

    ctr, fnctr = 0, 0
    
    # start with npts initial intervals    
    if no_pts > 1
        ints = Interval[]
        d = (b-a)/no_pts
        l = float(a)
        fl = f(l)
        fnctr += 1
        for i in 1:(no_pts-1)
            r = a + i * d
            m = midpoint(l, r)
            fm, fr = f(m), f(r)
            fnctr += 2
            push!(ints, Interval(l, m, r , fl, fm, fr))
            l,fl = r, fr
        end
        r = float(b)
        m = midpoint(l, r)
        fm, fr = f(m), f(r)
        fnctr += 2
        push!(ints, Interval(l, m, r, fl, fm, fr))
    else
        ints = Interval[Interval(f, float(a), float(b))]
        fnctr += 3
    end

    # container for zeros
    T = eltype(float(a))
    xzeros = T[]

    while !isempty(ints)
        ctr += 1
        alpha = pop!(ints)
        dx = deltax(alpha)

        # is the interval short? If so, find the root or discard
        lalpha, ralpha = split_interval(f, alpha)
        fnctr +=2

        if isbracket(alpha)
            # a bracket and small interval -> solve
            if dx <= BT(alpha, xatol, xrtol)
                push!(xzeros, find_zero(f, (alpha.l, alpha.r), Bisection()))
            continue
            end

            # keep the new bracket, check if other has promise
            if isbracket(lalpha) 
                push!(ints, lalpha)                    
                has_promise(ralpha, C, maxmultiplicity) && push!(ints, ralpha)
            else  
                push!(ints, ralpha)
                has_promise(lalpha, C, maxmultiplicity) && push!(ints, lalpha)
            end
            
        else
            
            # do we try to find zero in non-bracketing interval?
            # do we divide and try again?
            if isapproxzero(alpha, atol, rtol)
                # polish off to see if can find zero in interval
                try
                    rt = find_zero(f, alpha.m, Order2(), atol=atol, rtol=rtol)
                    if alpha.l < rt < alpha.r
                        push!(xzeros, rt)
                    end
                catch err
#                    push!(xzeros, alpha.m) # push anyways?
                end
                continue
            end

            ## this could check if lalpha or ralpha have promise
            ## has_promise(lalpha, C, maxmultiplicity) && push!(ints, lalpha)
            ## has_promise(rlapha, C, maxmultiplicity) && push!(ints, ralpha)            
            has_promise(alpha, lalpha, ralpha, C, maxmultiplicity) && append!(ints, (lalpha, ralpha))
            
        end

        # taking too long?
        ctr >= maxevals     && throw(ConvergenceFailed("too many evaluations $ctr $fnctr"))
        fnctr >= maxfnevals && throw(ConvergenceFailed("too many fn evaluations"))
        
    end
    
    if verbose
        nrts = length(xzeros)
        plural = nrts == 1 ? "" : "s"
        println("find_zeros found $nrts zero$plural in ($a, $b).")
        println(" * there were $ctr subintervals considered")
        println(" * there were $fnctr function evaluations taken")
        println("")
    end

    return xzeros
end
