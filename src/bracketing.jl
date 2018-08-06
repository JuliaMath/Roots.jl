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
function bisection64(f, a::T, b::T) where {T <: Union{Float64, Float32, Float16}}
    u, v = promote(float(a), float(b))
    if v < u
        u,v = v,u
    end

    isinf(u) && (u = nextfloat(u))
    isinf(v) && (u = prevfloat(u))

    su, sv = sign(u), sign(v)
    
    if su * sv < 0
        # move to 0
        c = zero(u)
        sfu, sfc = sign(f(u)), sign(f(c))
        if sfu == sfc
            u =  c
        else
            v = c
        end
    end
    
    T == Float64 && return _bisection64(T, UInt64, f, u, v)
    T == Float32 && return _bisection64(T, UInt32, f, u, v)
    return _bisection64(T, UInt16, f, u, v)
end

## a,b same sign or zero, sfa * sfb < 0 is assumed
function _bisection64(T, S, f, a, b)
    nan = (0*a)/(0*a)
    negate = sign(a) < 0 ? -one(T) : one(T)

    ai, bi = reinterpret(S, abs(a)), reinterpret(S, abs(b))
    
    fa = f(a)
    iszero(fa) && return a
    sfa = sign(f(a))
    iszero(f(b)) && return b
    ai == bi && return nan
    
    mi = (ai + bi ) >> 1
    m = negate * reinterpret(T, mi)

    while a < m < b
        

        sfm = sign(f(m))
        iszero(sfm) && return m
        isnan(sfm) && return m

        if sfa * sfm < 0
            b, bi = m,  mi
        else
            a, sfa, ai = m, sfm, mi
        end

        mi = (ai + bi) >> 1
        m = negate * reinterpret(T, mi)
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
there. This method avoids floating point issues and when the
tolerances are set to zero (the default) guarantees a "best" solution
(one where a zero is found or the bracketing interval is of the type
`[a, nextfloat(a)]`).

When tolerances are given, this algorithm terminates when the midpoint
is approximately equal to an endpoint using absolute tolerance `xatol`
and relative tolerance `xrtol`. 

When a zero tolerance is given and the values are not `Float64`
values, this will call the `A42` method which has guaranteed convergence.
    
"""
struct Bisection <: AbstractBisection end  # either solvable or A42
struct BisectionExact <: AbstractBisection end

"""
    Roots.A42()

Bracketing method which finds the root of a continuous function within
a provided interval [a, b], without requiring derivatives. It is based
on algorithm 4.2 described in: 1. G. E. Alefeld, F. A. Potra, and
Y. Shi, "Algorithm 748: enclosing zeros of continuous functions," ACM
Trans. Math. Softw. 21, 327–344 (1995), DOI: 10.1145/210089.210111 .
"""
mutable struct A42 <: AbstractBisection end

## tracks for bisection, different, we show bracketing interval
function log_step(l::Tracks, M::AbstractBisection, state)
    push!(l.xs, state.xn0)
    push!(l.xs, state.xn1) # we store [ai,bi, ai+1, bi+1, ...]
end
function show_tracks(l::Tracks, M::AbstractBisection)
    xs = l.xs
    n = length(xs)
    for (i,j) in enumerate(1:2:(n-1))
        println(@sprintf("(%s, %s) = (% 18.16f, % 18.16f)", "a_$(i-1)", "b_$(i-1)", xs[j], xs[j+1]))
    end
    println("")
end
        
    

## helper function
function adjust_bracket(x0)
    u, v = float.(promote(x0...))
    if u > v
        u,v = v,u
    end
    isinf(u) && (u = nextfloat(u))
    isinf(v) && (v = prevfloat(v))
    u, v
end

function init_state(method::AbstractBisection, fs, x)
    length(x) > 1 || throw(ArgumentError(bracketing_error))

    x0, x1 = adjust_bracket(x)
    y0, y1 = sign.(promote(fs(x0), fs(x1)))
    y0 * y1 > 0 && throw(ArgumentError("bracketing_error"))
    m = _middle(x0, x1)
    
    state = UnivariateZeroState(x1, x0, m,
                                y1, y0,
                                0, 2,
                                false, false, false, false,
                                "")
    state

end

function init_state!(state::UnivariateZeroState{T,S}, ::AbstractBisection, fs, x::Union{Tuple, Vector}) where {T, S}
    x0, x1 = adjust_bracket(x)
    fx0::S, fx1::S = sign(fs(x0)), sign(fs(x1))
    m = _middle(x0, x1)
    init_state!(state, x1, x0, m, fx1, fx0)
end

# for Bisection, the defaults are zero tolerances and strict=true
function init_options(::M,
                      state::UnivariateZeroState{T,S};
                      xatol=missing,
                      xrtol=missing,
                      atol=missing,
                      rtol=missing,
                      maxevals::Int=typemax(Int),
                      maxfnevals::Int=typemax(Int)) where {M <: Union{Bisection, BisectionExact,  A42}, T, S}

    ## Where we set defaults
    x1 = real(oneunit(state.xn1))
    fx1 = real(oneunit(float(state.fxn1)))
    strict = true

    # all are 0 by default
    options = UnivariateZeroOptions(ismissing(xatol) ? zero(x1) : xatol,       # unit of x
                                    ismissing(xrtol) ? zero(x1/oneunit(x1)) : xrtol,               # unitless
                                    ismissing(atol)  ? zero(fx1) : atol,  # units of f(x)
                                    ismissing(rtol)  ? zero(fx1/oneunit(fx1)) : rtol,            # unitless
                                    maxevals, maxfnevals, strict)

    options
end

function init_options!(options::UnivariateZeroOptions{Q,R,S,T}, ::Bisection) where {Q, R, S, T}
    options.xabstol = zero(Q)
    options.xreltol = zero(R)
    options.abstol = zero(S)
    options.reltol = zero(T)
    options.maxevals = typemax(Int)
    options.strict = true
end

## This uses _middle bisection Find zero using modified bisection
## method for FloatXX arguments.  This is guaranteed to take no more
## steps the bits of the type. The `a42` alternative usually has fewer
## iterations, but this seems to find the value with fewer function
## evaluations.
##
## This terminates when there is no more subdivision or function is zero

_middle(x::Float64, y::Float64) = _middle(Float64, UInt64, x, y)
_middle(x::Float32, y::Float32) = _middle(Float32, UInt32, x, y)
_middle(x::Float16, y::Float16) = _middle(Float16, UInt16, x, y)
_middle(x::Number, y::Number) = x + 0.5 * (y-x) # fall back or non Floats

function _middle(T, S, x, y)
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
    xint = reinterpret(S, abs(x))
    yint = reinterpret(S, abs(y))
    mid = (xint + yint) >> 1

    # reinterpret in original floating point
    unsigned = reinterpret(T, mid)

    negate ? -unsigned : unsigned
end

function update_state(method::Union{Bisection,BisectionExact}, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T<:Number,S<:Number}


    y0 = o.fxn0
    m::T = o.m  
    ym::S = sign(fs(m))
    incfn(o)

    if iszero(ym)
        o.message = "Exact zero found"
        o.xn1 = m 
        o.x_converged = true
        return nothing
    end
    
    if y0 * ym < 0
        o.xn1, o.fxn1 = m, ym
    else
        o.xn0, o.fxn0 = m, ym
    end

    o.m = _middle(o.xn0, o.xn1)            
    return nothing

end

## convergence is much different here
## the method converges,
## as we bound between x0, nextfloat(x0) is not measured by eps(), but eps(x0)
function assess_convergence(method::Union{Bisection}, state::UnivariateZeroState{T,S}, options) where {T, S}

   
    state.x_converged && return true

    x0, x1, m::T = state.xn0, state.xn1, state.m

    if !(x0 < m < x1)
        state.x_converged = true
        return true
    end

    tol = max(options.xabstol, max(abs(x0), abs(x1)) * options.xreltol)
    if x1 - x0 > tol 
        return false
    end
    
    
    state.message = ""
    state.x_converged = true
    return true
end

# for exact convergence, we can skip some steps
function assess_convergence(method::BisectionExact, state::UnivariateZeroState{T,S}, options) where {T, S}

    state.x_converged && return true
    
    x0, m::T, x1 = state.xn0, state.m, state.xn1

    x0 < m < x1 && return false

    state.x_converged = true
    return true
end



## Bisection has special cases
## for FloatNN types, we have a slightly faster `bisection64` method
## for zero tolerance, we have either BisectionExact or A42 methods
## for non-zero tolerances, we have either a general Bisection or an A42
function find_zero(fs, x0, method::M;
                   tracks = NullTracks(),
                   verbose=false,
                   kwargs...) where {M <: Union{Bisection, A42}}
    
    x = adjust_bracket(x0)
    T = eltype(x[1])
    F = callable_function(fs)
    state = init_state(method, F, x)
    options = init_options(method, state; kwargs...)
    tol = max(options.xabstol, maximum(abs.(x)) * options.xreltol)

    l = (verbose && isa(tracks, NullTracks)) ? Tracks(eltype(state.xn1)[], eltype(state.fxn1)[]) : tracks
    
    if iszero(tol)
        if T <: FloatNN
            !verbose && return bisection64(F, state.xn0, state.xn1) # speedier
            find_zero(BisectionExact(), F, options, state, l)
        else
            return a42(F, state.xn0, state.xn1, xtol=zero(T), verbose=verbose)
        end
    else
        find_zero(method, F, options, state, l)
    end

    if verbose
        show_trace(method, state, l)
    end
    
    state.xn1
    
end

## The Roots.A42() method is not implemented within the frame work
function find_zero(method::A42, F, options::UnivariateZeroOptions, state::UnivariateZeroState{T,S}, tracks) where {T<:Number, S<:Number}
    x0, x1 = state.xn0, state.xn1
    tol = max(options.xabstol, max(abs(x0), abs(x1)) * options.xreltol)
    state.xn1 = a42(F, x0, x1; xtol=tol, maxeval=options.maxevals),
    state.message = "Used Alefeld-Potra-Shi method, `Roots.a42`, to find the zero. Iterations and function evaluations are not counted properly."
    state.stopped = state.x_converged  = true
    
    return state.xn1
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

    u, v = adjust_bracket((a,b))
    fu, fv = f(u), f(v)

    if sign(fu)*sign(fv) >= 0
        error("on input a < b and f(a)f(b) < 0 must both hold")
    end
    
    if xtol/oneunit(xtol) < 0.0
        error("tolerance must be >= 0.0")
    end
    
    c, fc = secant_step(f, u, fu, v, fv)
    a42a(f, u, fu, v, fv, float(c), fc,
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
    2u*eps(one(T)) + tol
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
    if iszero(fc)
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
function secant_step(f, a::T, fa, b, fb) where {T}

    c = a - fa/(fb - fa)*(b - a)
    tol = 5*eps(one(a))
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
    if iszero(A)
        return secant_step(f, a, fa, b, fb)
    end
    
    r = A*fa/oneunit(A*fa) > 0 ? a : b
    for i = 1:k
        r -= (fa + (B + A*(r - b))*(r - a))/(B + A*(2*r - a - b))
    end
    if isnan(r) || (r <= a || r >= b)
        r, fr = secant_step(f, a, fa, b, fb)
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
struct FalsePosition{R} <: AbstractBisection end
FalsePosition(x=:anderson_bjork) = FalsePosition{x}()

function update_state(method::FalsePosition, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}

    a::T, b::T =  o.xn0, o.xn1

    fa::S, fb::S = o.fxn0, o.fxn1

    lambda = fb / (fb - fa)
    tau = 1e-10                   # some engineering to avoid short moves
    if !(tau < abs(lambda) < 1-tau)
        lambda = 1/2
    end
    x::T = b - lambda * (b-a)        
    fx::S = fs(x)
    incfn(o)

    if iszero(fx)
        o.xn1 = x
        o.fxn1 = fx
        return
    end

    if sign(fx)*sign(fb) < 0
        a, fa = b, fb
    else
        fa = galdino_reduction(method, fa, fb, fx) #galdino[method.reduction_factor](fa, fb, fx)
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

# from Chris Elrod; https://raw.githubusercontent.com/chriselrod/AsymptoticPosteriors.jl/master/src/false_position.jl
@generated function galdino_reduction(methods::FalsePosition{R}, fa, fb, fx) where {R}
    f = galdino[R]
    quote
        $Expr(:meta, :inline)
        $f(fa, fb, fx)
    end
end
