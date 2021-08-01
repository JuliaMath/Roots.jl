###


const bracketing_error = """The interval [a,b] is not a bracketing interval.
You need f(a) and f(b) to have different signs (f(a) * f(b) < 0).
Consider a different bracket or try fzero(f, c) with an initial guess c.

"""



abstract type AbstractBisection <: AbstractBracketing end
fn_argout(::AbstractBracketing) = 1

"""

    Bisection()
    Roots.BisectionExact()

If possible, will use the bisection method over `Float64` values. The
bisection method starts with a bracketing interval `[a,b]` and splits
it into two intervals `[a,c]` and `[c,b]`, If `c` is not a zero, then
one of these two will be a bracketing interval and the process
continues. The computation of `c` is done by `_middle`, which
reinterprets floating point values as unsigned integers and splits
there. It was contributed  by  [Jason Merrill](https://gist.github.com/jwmerrill/9012954).
This method avoids floating point issues and when the
tolerances are set to zero (the default) guarantees a "best" solution
(one where a zero is found or the bracketing interval is of the type
`[a, nextfloat(a)]`).

When tolerances are given, this algorithm terminates when the midpoint
is approximately equal to an endpoint using absolute tolerance `xatol`
and relative tolerance `xrtol`.

When a zero tolerance is given and the values are not `Float64`
values, this will call the [`A42`](@ref) method.


"""
struct Bisection <: AbstractBisection end  # either solvable or A42
struct BisectionExact <: AbstractBisection end

## tracks for bisection, different, we show bracketing interval
function log_step(l::Tracks, M::AbstractBracketing, state)
    push!(l.xs, state.xn0)
    push!(l.xs, state.xn1) # we store [ai,bi, ai+1, bi+1, ...]
end
function show_tracks(l::Tracks, M::AbstractBracketing)

    xs = l.xs
    n = length(xs)
    for (i,j) in enumerate(1:2:(n-1))
        println(@sprintf("(%s, %s) = (% 18.16f, % 18.16f)", "a_$(i-1)", "b_$(i-1)", xs[j], xs[j+1]))
    end
    println("")
end



## helper function
function adjust_bracket(x0)
    u, v = float.(promote(_extrema(x0)...))
    if u > v
        u,v = v,u
    end
    isinf(u) && (u = nextfloat(u))
    isinf(v) && (v = prevfloat(v))
    u, v
end


function init_state(M::AbstractBracketing, F::Callable_Function, x; kwargs...)
    x₀, x₁ = adjust_bracket(x)
    fx₀, fx₁ = promote(float(first(F(x₀))), float(first(F(x₁))))
    state::UnivariateZeroState = init_state(M, x₀, x₁, fx₀, fx₁;
                                            fnevals = initial_fnevals(M),
                                            kwargs...)
    init_state!(state, M, F, clear=false)
    state
end
initial_fnevals(::Roots.AbstractBracketing) = 2

function init_state!(state, M::AbstractBisection, F::Callable_Function; clear=true)
    x0, x1 = state.xn0, state.xn1
    fx0, fx1 = state.fxn0, state.fxn1

    if x0 > x1
        x0, x1, fx0, fx1 = x1, x0, fx1, fx0
    end

    fx0fx1 = sign(fx0) * sign(fx1)
    if iszero(fx0fx1)
        # should be an error -- not bracketing, but we have a zero, so return it.
        m,fm = iszero(fx0) ? (x0,fx0) : (x1, fx1)
        state.x_converged = true
        state.f_converged = true
        state.xstar  = m
        state.fxstar =  fm
        return state
    end

    if isnan(fx0fx1)
        m,fm = isnan(fx0) ? (x0,fx0) : (x1, fx1)
        state.x_converged = true
        state.f_converged = true
        state.message = "NaN encountered. "
        state.xstar  = m
        state.fxstar =  fm
        return state
    end

    assert_bracket(fx0, fx1)

    # we want a,b to be same sign, finite
    if sign(x0) * sign(x1) < 0
        m = zero(x1)
        fm = F(m) # iszero(fm) caught in assess_convergence
        incfn(state)
        if sign(fx0) * sign(fm) < 0
            x1, fx1 = m, fm
        else  sign(fx0) * sign(fm) > 0
            x0, fx0 = m, fm
        end
    end

    m = __middle(x0, x1)
    fm = F(m)
    incfn(state)

    state.xn0, state.xn1 = x0, x1
    state.fxn0, state.fxn1 = fx0, fx1
    empty!(state.m); empty!(state.fm)
    push!(state.m, m), push!(state.fm, fm)

    clear && clear_convergence_flags!(state)
    nothing
end



# for Bisection, the defaults are zero tolerances and strict=true
"""
    default_tolerances(M::Bisection, [T], [S])


For `Bisection` (or `BisectionExact`), when the `x` values are of type `Float64`, `Float32`,
or `Float16`, the default tolerances are zero and there is no limit on
the number of iterations or function evalutions. In this case, the
algorithm is guaranteed to converge to an exact zero, or a point where
the function changes sign at one of the answer's adjacent floating
point values.

For other types,  the [`A42`](@ref) method (with its tolerances) is used.

"""
default_tolerances(M::Union{Bisection, BisectionExact}) = default_tolerances(M,Float64, Float64)
function default_tolerances(::M, ::Type{T}, ::Type{S}) where {M<:Union{Bisection, BisectionExact},T,S}
    xatol = zero(T)
    xrtol = zero(one(T))
    atol = zero(float(one(S))) * oneunit(S)
    rtol = zero(float(one(S))) * one(S)
    maxevals = typemax(Int)
    maxfnevals = typemax(Int)
    strict = true
    (xatol, xrtol, atol, rtol, maxevals, maxfnevals, strict)
end


# find middle of (a,b) with convention that
# * if a, b finite, they are made non-finite
# if a,b of different signs, middle is 0
# middle falls back to a/2 + b/2, but
# for Float64 values, middle is over the
# reinterpreted unsigned integer.
function _middle(x, y)
    a = isinf(x) ? nextfloat(x) : x
    b = isinf(y) ? prevfloat(y) : y
    if sign(a) * sign(b) < 0
        return zero(a)
    else
        __middle(a, b)
    end
end

const FloatNN = Union{Float64, Float32, Float16}

## find middle assuming a,b same sign, finite
## Alternative "mean" definition that operates on the binary representation
## of a float. Using this definition, bisection will never take more than
## 64 steps (over Float64)
__middle(x::Float64, y::Float64) = __middle(Float64, UInt64, x, y)
__middle(x::Float32, y::Float32) = __middle(Float32, UInt32, x, y)
__middle(x::Float16, y::Float16) = __middle(Float16, UInt16, x, y)
## fallback for non FloatNN number types
__middle(x::Number, y::Number) = 0.5*x + 0.5*y


function __middle(T, S, x, y)
    # Use the usual float rules for combining non-finite numbers
    # do division over unsigned integers with bit shift
    xint = reinterpret(S, abs(x))
    yint = reinterpret(S, abs(y))
    mid = (xint + yint) >> 1

    # reinterpret in original floating point
    sign(x+y) * reinterpret(T, mid)
end




## the method converges,
## as we bound between x0, nextfloat(x0) is not measured by eps(), but eps(x0)
function assess_convergence(M::Bisection, state::UnivariateZeroState{T,S}, options) where {T, S}

    assess_convergence(BisectionExact(), state, options) && return true
    x0, x1, xm = state.xn0, state.xn1, first(state.m)

    xtol = max(options.xabstol, max(abs(x0), abs(x1)) * options.xreltol)
    x_converged = x1 - x0 < xtol

    if x_converged
        state.message=""
        state.xstar = xm
        state.x_converged = true
        return true
    end


    fm = first(state.fm)
    ftol = max(options.abstol, max(abs(x0), abs(x1)) * options.reltol)
    f_converged = abs(fm) < ftol

    if f_converged
        state.message = ""
        state.xstar = xm
        state.fxstar = fm
        state.f_converged = f_converged
        return true
    end

    return false
end

# for exact convergence, we can skip some steps
function assess_convergence(M::BisectionExact, state::UnivariateZeroState{T,S}, options) where {T, S}

    state.x_converged && return true

    x0, xm::T, x1 = state.xn0, state.m[1], state.xn1
    y0, ym, y1 = state.fxn0, state.fm[1], state.fxn1

    for (c,fc) in ((x0,y0), (xm,ym), (x1, y1))
        if iszero(fc) || isnan(fc) #|| isinf(fc)
            state.f_converged = true
            state.x_converged = true
            state.xstar = c
            state.fxstar = fc
            return true
        end
    end

    x0 < xm < x1 && return false

    state.xstar = x1
    state.fxstar = ym
    state.x_converged = true
    return true
end



##################################################


function update_state(method::Union{Bisection, BisectionExact}, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T<:Number,S<:Number}

    y0 = o.fxn0
    m::T = o.m[1]
    ym::S = o.fm[1] #fs(m)

    if sign(y0) * sign(ym) < 0
        o.xn1, o.fxn1 = m, ym
    else
        o.xn0, o.fxn0 = m, ym
    end

    m  = __middle(o.xn0, o.xn1) # assume a,b have same sign
    fm::S = fs(m)
    o.m[1], o.fm[1] = m, fm
    incfn(o)

    return nothing

end

##################################################

## Bisection has special cases
## for zero tolerance, we have either BisectionExact or A42 methods
## for non-zero tolerances, we have use thegeneral Bisection method
function find_zero(fs, x0, method::M;
                   p=nothing,
                   tracks = NullTracks(),
                   verbose=false,
                   kwargs...) where {M <: Union{Bisection}}

    F = Callable_Function(method, fs, p)
    state = init_state(method, F, x0)
    options = init_options(method, state; kwargs...)
    l = Tracks(verbose, tracks, state)

    # check if tolerances are exactly 0
    iszero_tol = iszero(options.xabstol) && iszero(options.xreltol) && iszero(options.abstol) && iszero(options.reltol)

    # zero_tol and use either exact or A42
    iszero_tol && (method = xType(state) <: FloatNN ? BisectionExact() : A42())

    ZPI = ZeroProblemIterator(method, F, state, options, l)
    solve!(ZPI)

    verbose && tracks(ZPI)

    state.xstar

end

## --------------------------------------------------
##
# assume int[xn0,xn1] is a bracket in state
function assess_convergence(method::AbstractBracketing, state::UnivariateZeroState{T,S}, options) where {T,S}

    if state.stopped || state.x_converged || state.f_converged
        return true
    end

    check_steps_fnevals(state, options) && return true

    a,b,fa,fb = state.xn0, state.xn1, state.fxn0, state.fxn1

    if isnan(a) || isnan(b)
        state.convergence_failed = true
        state.message *= "NaN produced by algorithm. "
        return true
    end

    M = maximum(abs, (a,b))
    δₓ = maximum(promote(options.xabstol, M * options.xreltol, sign(options.xreltol) *   eps(M)))

    if abs(b-a) <= 2δₓ
        state.x_converged = true
    end

    # check f
    u,fu = choose_smallest(a,b,fa,fb)
    δ = maximum(promote(options.abstol, M * options.reltol * (oneunit(fu) / oneunit(u))))
    if abs(fu) <= δ
        state.f_converged = true
        if iszero(fu)
            state.x_converged
            state.message *= "Exact zero found. "
        end
    end

    if state.f_converged || state.x_converged
        state.xstar = u
        state.fxstar = fu
        return true
    end

    return false
end

## convergence is much different here
function check_zero(::AbstractBracketing, state, c, fc)
    if isnan(c)
        state.stopped = true
        #state.xn1 = c
        state.xstar = c
        state.fxstar =fc
        state.message *= "NaN encountered. "
        return true
    elseif isinf(c)
        state.stopped = true
        #state.xn1 = c
        state.xstar = c
        state.fxstar =fc
        state.message *= "Inf encountered. "
        return true
    elseif iszero(fc)
        state.f_converged=true
        state.message *= "Exact zero found. "
        state.xstar = c
        state.fxstar =  fc
        return true
    end
    return false
end

###################################################
#
## Alefeld, Potra, Shi have two algorithms belosw, one is most efficient, but
## slightly slower than other.
abstract type AbstractAlefeldPotraShi <: AbstractBracketing end

"""
    Roots.A42()

Bracketing method which finds the root of a continuous function within
a provided interval [a, b], without requiring derivatives. It is based
on algorithm 4.2 described in: G. E. Alefeld, F. A. Potra, and
Y. Shi, "Algorithm 748: enclosing zeros of continuous functions," ACM
Trans. Math. Softw. 21, 327–344 (1995), DOI: [10.1145/210089.210111](https://doi.org/10.1145/210089.210111).
Originally by John Travers.

"""
struct A42 <: AbstractAlefeldPotraShi end

## put in utils?
@inline isbracket(fa,fb) = sign(fa) * sign(fb) < 0
assert_bracket(fx0, fx1) = isbracket(fx0, fx1) || throw(ArgumentError(bracketing_error))

# f[a, b]
@inline f_ab(a,b,fa,fb) = (fb - fa) / (b-a)

# f[a,b,d]
@inline function f_abd(a,b,d,fa,fb,fd)
    fab, fbd = f_ab(a,b,fa,fb), f_ab(b,d,fb,fd)
    (fbd - fab)/(d-a)
end

# a bit better than a - fa/f_ab
@inline secant_step(a, b, fa, fb) =  a - fa * (b - a) / (fb - fa)

# assume fc != 0
## return a1,b1,d with a < a1 <  < b1 < b, d not there
@inline function bracket(a,b,c, fa, fb, fc)
    if isbracket(fa, fc)
        # switch b,c
        return (a,c,b, fa, fc, fb)
    else
        # switch a,c
        return (c,b,a, fc, fb, fa)
    end
end

# Cubic if possible, if not, quadratic(3)
function take_a42_step(a::T, b, d, ee, fa, fb, fd, fe, k, delta=zero(T)) where {T}

    # if r is NaN or Inf we move on by condition. Faster than checking ahead of time for
    # distinctness
    r = ipzero(a,b,d,ee, fa, fb,fd,fe, delta) # let error and see difference in allcoation?
    (a + 2delta < r < b - 2delta) && return r
    r = newton_quadratic(a,b,d,fa,fb,fd, 3, delta)
end

function ipzero(a::T, b, c, d, fa, fb, fc, fd, delta=zero(T)) where {T}
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

    (a + 2delta < c < b - 2delta) && return c

    newton_quadratic(a,b,d,fa,fb,fd, 3, delta)

end

# return c in (a+delta, b-delta)
# adds part of `bracket` from paper with `delta`
function newton_quadratic(a::T, b, d, fa, fb, fd, k::Int, delta=zero(T)) where {T}

    A = f_abd(a,b,d,fa,fb,fd)
    r = isbracket(A,fa) ? b : a

    # use quadratic step; if that fails, use secant step; if that fails, bisection
    if !(isnan(A) || isinf(A)) || !iszero(A)
        B = f_ab(a,b,fa,fb)

        for i in 1:k
            Pr = fa + B * (r-a) +  A * (r-a)*(r-b)
            Prp = (B + A*(2r - a - b))
            r -= Pr / Prp
        end
        if a+2delta < r < b - 2delta
            return r
        end
    end
    # try secant step
    r =  secant_step(a, b, fa, fb)

    if a + 2delta < r < b - 2delta
        return r
    end

    return _middle(a, b) # is in paper r + sgn * 2 * delta

end

function init_state!(state::UnivariateZeroState{T,S}, M::AbstractAlefeldPotraShi, F::Callable_Function;
                     middle=nothing, clear=true) where {T,S}


    a,b,fa,fb = state.xn0, state.xn1, state.fxn0, state.fxn1

    if a > b
        a,b,fa,fb = b,a,fb,fa
    end



    c::T = (middle != nothing && a < middle < b) ? middle : _middle(a,b)
    fc::S = F(c)
    incfn(state)

    a,b,d,fa,fb,fd = bracket(a,b,c,fa,fb,fc)
    ee, fe = d, fd

    state.xn0, state.xn1, state.fxn0, state.fxn1 = a,b,fa,fb
    empty!(state.m), empty!(state.fm)
    append!(state.m, (d,ee)), append!(state.fm, (fd, fe))

    clear && clear_convergence_flags!(state)

    nothing
end



# for A42, the defaults are reltol=eps(), atol=0; 45 evals and strict=true
# this *basically* follows the tol in the paper (2|u|*rtol + atol)
"""
    default_tolerances(::AbstractAlefeldPotraShi, T, S)

The default tolerances for Alefeld, Potra, and Shi methods are
`xatol=zero(T)`, `xrtol=eps(T)/2`, `atol= zero(S), and rtol=zero(S)`, with
appropriate units; `maxevals=45`, `maxfnevals = Inf`; and `strict=true`.

"""
default_tolerances(M::AbstractAlefeldPotraShi) = default_tolerances(M, Float64, Float64)
function default_tolerances(::AbstractAlefeldPotraShi, ::Type{T}, ::Type{S}) where {T,S}
    xatol = zero(T)
    xrtol = eps(one(T))/2
    atol = zero(float(one(S))) * oneunit(S)
    rtol = zero(float(one(S))) * one(S)
    maxevals = 45
    maxfnevals = typemax(Int)
    strict = true
    (xatol, xrtol, atol, rtol, maxevals, maxfnevals, strict)
end





## initial step, needs to log a,b,d
function log_step(l::Tracks, M::AbstractAlefeldPotraShi, state, ::Any)
    a, b, c = state.xn0, state.xn1, state.m[1]
    append!(l.xs, extrema((a,b,c)))
    push!(l.xs, a)
    push!(l.xs, b) # we store [ai,bi, ai+1, bi+1, ...] for brackecting methods
end

# Main algorithm for A42 method
function update_state(M::A42, f, state::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}

    a::T,b::T,d::T, ee::T = state.xn0, state.xn1, state.m[1], state.m[2]
    fa::S,fb::S,fd::S,fe::S = state.fxn0, state.fxn1, state.fm[1], state.fm[2]

    an, bn = a, b
    μ, λ = 0.5, 0.7
    tole = max(options.xabstol, max(abs(a),abs(b)) * options.xreltol) # paper uses 2|u|*rtol + atol
    delta = λ * tole

    if state.steps < 1
        c = newton_quadratic(a, b, d, fa, fb, fd, 2)
    else
        c = ipzero(a, b, d, ee, fa, fb, fd, fe)
    end
    fc::S = f(c)
    incfn(state)
    if check_zero(M, state, c, fc)
        return nothing
    end

    ab::T, bb::T, db::T, fab::S, fbb::S, fdb::S = bracket(a,b,c,fa,fb,fc)
    eb::T, feb::S = d, fd

    cb::T = take_a42_step(ab, bb, db, eb, fab, fbb, fdb, feb, delta)
    fcb::S = f(cb)
    incfn(state)
    if check_zero(M, state, cb, fcb)
        # tighten up bracket
        state.xn0, state.xn1, state.m[1]  = ab, bb, db
        state.fxn0, state.fxn1, state.fm[1]= fab, fbb, fdb

        return nothing
    end

    ab,bb,db,fab,fbb,fdb = bracket(ab,bb,cb,fab,fbb,fcb)


    u::T, fu::S = choose_smallest(ab, bb, fab, fbb)
    cb = u - 2 * fu * (bb - ab) / (fbb - fab)
    ch::T = cb
    if abs(cb - u) > 0.5 * (b-a)
        ch = _middle(an, bn)
    end
    fch::S = f(ch)
    incfn(state)
    if check_zero(M, state, ch, fch)
        # tighten up bracket
        state.xn0, state.xn1, state.m[1]  = ab, bb, db
        state.fxn0, state.fxn1, state.fm[1]= fab, fbb, fdb

        return nothing
    end

    ah::T, bh::T, dh::T, fah::S, fbh::S, fdh::S = bracket(ab, bb, ch, fab, fbb, fch)

    if bh - ah < μ * (b - a)
        #a, b, d, fa, fb, fd = ahat, b, dhat, fahat, fb, fdhat # typo in paper
        a, b, d, ee =  ah, bh, dh, db
        fa, fb, fd, fe = fah, fbh, fdh, fdb
    else
        m::T = _middle(ah, bh)
        fm::S = f(m)
        incfn(state)
        ee, fe = dh, fdh
        a, b, d, fa, fb, fd = bracket(ah, bh, m, fah, fbh, fm)
    end
    state.xn0, state.xn1, state.m[1], state.m[2],  = a, b, d, ee
    state.fxn0, state.fxn1, state.fm[1], state.fm[2] = fa, fb, fd, fe

    return nothing
end


####
"""
    Roots.AlefeldPotraShi()

Follows algorithm in "ON ENCLOSING SIMPLE ROOTS OF NONLINEAR
EQUATIONS", by Alefeld, Potra, Shi; DOI:
[10.1090/S0025-5718-1993-1192965-2](https://doi.org/10.1090/S0025-5718-1993-1192965-2).
 Efficiency is 1.618. Less efficient, but can be faster than the [`A42`](@ref) method.
Originally by John Travers.
"""
struct AlefeldPotraShi <: AbstractAlefeldPotraShi end

# ## 3, maybe 4, functions calls per step
function update_state(M::AlefeldPotraShi, f, state::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}

    a::T,b::T,d::T = state.xn0, state.xn1, state.m[1]
    fa::S,fb::S,fd::S = state.fxn0, state.fxn1, state.fm[1]

    μ, λ = 0.5, 0.7
    tole = max(options.xabstol, max(abs(a),abs(b)) * options.xreltol) # paper uses 2|u|*rtol + atol
    delta = λ * tole

    c::T = newton_quadratic(a, b, d, fa, fb, fd, 2, delta)
    fc::S = f(c)
    incfn(state)
    check_zero(M, state, c, fc) && return nothing

    a,b,d,fa,fb,fd = bracket(a,b,c,fa,fb,fc)

    c = newton_quadratic(a,b,d,fa,fb,fd, 3, delta)
    fc = f(c)
    incfn(state)
    if check_zero(M, state, c, fc)
        # tighten up bracket
        state.xn0, state.xn1, state.m[1]  = a, b, d
        state.fxn0, state.fxn1, state.fm[1]= fa, fb, fd

        return nothing
    end

    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb,fc)

    u::T, fu::S = choose_smallest(a, b, fa, fb)
    c = u - 2 * fu * (b - a) / (fb - fa)
    if abs(c - u) > 0.5 * (b - a)
        c = _middle(a, b)
    end
    fc = f(c)
    incfn(state)
    if  check_zero(M, state, c, fc)
        # tighten up bracket
        state.xn0, state.xn1, state.m[1]  = a, b, d
        state.fxn0, state.fxn1, state.fm[1]= fa, fb, fd

        return nothing
    end

    ahat::T, bhat::T, dhat::T, fahat::S, fbhat::S, fdhat::S = bracket(a, b, c, fa, fb, fc)
    if bhat - ahat < μ * (b - a)
        #a, b, d, fa, fb, fd = ahat, b, dhat, fahat, fb, fdhat # typo in paper
        a, b, d, fa, fb, fd = ahat, bhat, dhat, fahat, fbhat, fdhat
    else
        m::T = _middle(ahat, bhat)
        fm::S = f(m)
        incfn(state)
        a, b, d, fa, fb, fd = bracket(ahat, bhat, m, fahat, fbhat, fm)
    end
    state.xn0, state.xn1, state.m[1] = a, b, d
    state.fxn0, state.fxn1, state.fm[1] = fa, fb, fd

    return nothing
end


### Brent
"""
    Roots.Brent()

An implementation of
[Brent's](https://en.wikipedia.org/wiki/Brent%27s_method) (or Brent-Dekker) method.
This method uses a choice of inverse quadratic interpolation or a secant
step, falling back on bisection if necessary.

"""
struct Brent <: AbstractBracketing end

function log_step(l::Tracks, M::Brent, state)
    a,b = state.xn0, state.xn1
    u,v = a < b ? (a,b) : (b,a)
    push!(l.xs, a)
    push!(l.xs, b) # we store [ai,bi, ai+1, bi+1, ...]
end


# we store mflag as -1, or +1 in state.m[2]
function init_state!(state::UnivariateZeroState{T,S}, ::Brent, F; clear=true) where {T,S}

    u, v, fu, fv = state.xn0, state.xn1, state.fxn0, state.fxn1

    empty!(state.m); empty!(state.fm)
    # brent store b as smaller of |fa|, |fb|
    if abs(fu) > abs(fv)
        state.xn0, state.xn1, state.fxn0, state.fxn1 = u, v, fu, fv
        append!(state.m, (u,u)); append!(state.fm, (fu, one(fu)))
    else
        state.xn0, state.xn1, state.fxn0, state.fxn1 = v, u, fv, fu
        append!(state.m, (v,v)); append!(state.fm, (fv, one(fv)))
    end
    clear && clear_convergence_flags!(state)
    return nothing
end


function update_state(M::Brent, f, state::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}
    mflag = state.fm[2] > 0
    a, b, c, d = state.xn0, state.xn1, state.m[1], state.m[2]
    fa, fb, fc = state.fxn0, state.fxn1, state.fm[1]

    # next setp
    s::T = zero(a)
    if  !iszero(fa-fc) && !iszero(fb-fc)
        s =  a * fb * fc / (fa - fb) / (fa - fc) # quad step
        s += b * fa * fc / (fb - fa) / (fb - fc)
        s += c * fa * fb / (fc - fa) / (fc - fb)
        s
    else
        s = secant_step(a,b,fa,fb)
    end
    fs::S = f(s)
    incfn(state)
    check_zero(M, state, s, fs) && return nothing

    # guard step
    u,v = (3a+b)/4, b
    if u > v
        u,v = v, u
    end

    tol = max(options.xabstol, max(abs(b), abs(c), abs(d)) * options.xreltol)
    if !(u < s < v) ||
        (mflag && abs(s - b) >= abs(b-c)/2) ||
        (!mflag && abs(s - b) >= abs(b-c)/2) ||
        (mflag && abs(b-c) <= tol) ||
        (!mflag && abs(c-d) <= tol)
        s = _middle(a, b)
        fs = f(s)
        incfn(state)
        check_zero(M, state, s, fs) && return nothing
        mflag = true
    else
        mflag = false
    end

    d = c
    c,fc = b,fb

    if sign(fa) * sign(fs) < 0
        b, fb = s, fs
    else
        a, fa = s, fs
    end

    if abs(fa) < abs(fb)
        a, b, fa, fb = b, a, fb, fa
    end

    state.xn0, state.xn1, state.m[1], state.m[2] = a, b, c, d
    state.fxn0, state.fxn1, state.fm[1] = fa, fb, fc
    state.fm[2] = mflag ? one(fa) : -one(fa)

    return nothing
end


## ----------------------------

struct FalsePosition{R} <: AbstractBisection end
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
find_zero(x -> x^5 - x - 1, (-2, 2), FalsePosition())
```
"""
FalsePosition
FalsePosition(x=:anderson_bjork) = FalsePosition{x}()

# use fallback for derivative free
function assess_convergence(method::FalsePosition, state::UnivariateZeroState{T,S}, options) where {T,S}
    assess_convergence(Any, state, options)
end

function update_state(method::FalsePosition, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}

    a::T, b::T =  o.xn0, o.xn1
    fa::S, fb::S = o.fxn0, o.fxn1

    lambda = fb / (fb - fa)
    tau = 1e-10                   # some engineering to avoid short moves; still fails on some
    if !(tau < abs(lambda) < 1-tau)
        lambda = 1/2
    end

    x::T = b - lambda * (b-a)
    fx::S = fs(x)
    incfn(o)
    if iszero(fx)
        o.xn1 = x
        o.fxn1 = fx
        o.f_converged = true
        return
    end

    if sign(fx)*sign(fb) < 0
        a, fa = b, fb
    else
        fa = galdino_reduction(method, fa, fb, fx)
    end
    b, fb = x, fx

    o.xn0, o.xn1 = a, b
    o.fxn0, o.fxn1 = fa, fb

    nothing
end

# the 12 reduction factors offered by Galadino
# In RootsTesting.jl, we can see :12 has many more failures.
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




"""
    find_bracket(f, x0, method=A42(); kwargs...)

For bracketing methods returns an approximate root, the last bracketing interval used, and a flag indicating if an exact zero was found as a named tuple.

With the default tolerances, one  of these should be the case: `exact`  is `true` (indicating termination  of the algorithm due to an exact zero  being identified) or the length of `bracket` is less or equal than `2eps(maximum(abs.(bracket)))`. In the `BisectionExact` case, the 2 could  be replaced by 1, as the bracket, `(a,b)` will satisfy  `nextfloat(a) == b `;  the Alefeld,  Potra, and Shi algorithms don't quite have  that promise.

"""
function find_bracket(fs, x0, method::M=A42(); kwargs...) where {M <: Union{AbstractAlefeldPotraShi, BisectionExact}}
    x = adjust_bracket(x0)
    F = Callable_Function(method, fs) #callable_function(fs)
    state = init_state(method, F, x)
    options = init_options(method, state; kwargs...)

    # check if tolerances are exactly 0
    iszero_tol = iszero(options.xabstol) && iszero(options.xreltol) && iszero(options.abstol) && iszero(options.reltol)

    ZPI = ZeroProblemIterator(method, F, state, options, NullTracks())
    solve!(ZPI)

    (xstar=state.xstar, bracket=(state.xn0, state.xn1), exact=iszero(state.fxstar))

end
