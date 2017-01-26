###

# type to throw on succesful convergence
type StateConverged
    x0::Real
end

# type to throw on failure
type ConvergenceFailed
    reason::AbstractString
end

bracketing_error = """The interval [a,b] is not a bracketing interval.
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
# 64 steps.

function _middle(x::Float64, y::Float64)
  # Use the usual float rules for combining non-finite numbers
  if !isfinite(x) || !isfinite(y)
    return x + y
  end

  # Always return 0.0 when inputs have opposite sign
  if sign(x) != sign(y) && x != 0.0 && y != 0.0
    return 0.0
  end

  negate = x < 0.0 || y < 0.0

  xint = reinterpret(UInt64, abs(x))
  yint = reinterpret(UInt64, abs(y))
  unsigned = reinterpret(Float64, (xint + yint) >> 1)

  return negate ? -unsigned : unsigned
end

## fall back or non Floats
function _middle(x::Real, y::Real)
    x + (y-x)/2
end


####
## find_zero interface. We need to specialize for T<:Float64, and BigSomething
const BigSomething = Union{BigFloat, BigInt}

@compat abstract type AbstractBisection <: UnivariateZeroMethod end
type Bisection <: AbstractBisection end
type A42 <: AbstractBisection end

## This is a bit clunky, as we use `a42` for bisection when we don't have `Float64` values.
## As we don't have the `A42` algorithm implemented through `find_zero`, we adjust here.
function find_zero{M<:AbstractBisection, T<:Real}(f, x0::Vector{T}, method::M; maxevals::Int=50, verbose::Bool=false, kwargs...)
    x = sort(float(x0))
    if eltype(x) <: Float64
        prob, options = derivative_free_setup(method, DerivativeFree(f), x; verbose=verbose, maxevals=maxevals, kwargs...)
        find_zero(prob, method, options)
    else
        a42(f, x[1], x[2]; maxeval=maxevals, verbose=verbose)
    end
end

# call a42 in this case
function find_zero{T<:BigSomething, S}(method::Bisection, fs::DerivativeFree, o::UnivariateZeroState{T, S}, options::UnivariateZeroOptions)
    xn0, xn1 = sort([o.xn0, o.xn1])
    o.xn1 = a42(fs.f, o.xn0, o.xn1; xtol=options.xreltol, maxeval=options.maxevals, verbose=options.verbose)
    o.message = "Used Alefeld-Potra-Shi method, `Roots.a42`, to find the zero. Iterations and function evaluations are not counted properly."
    o.x_converged = true

    nothing
end

## This uses _middle bisection
## Find zero using modified bisection method for Float64 arguments.
## This is guaranteed to take no more than 64 steps. The `a42` alternative usually has
## fewer iterations, but this seems to find the value with fewer function evaluations.
##
## Terminates with `x1` when the bracket length of `[x0,x2]` is `<= max(xtol, xtolrel*abs(x1))` where `x1` is the midpoint .
## The tolerances can be set to 0, in which case, the termination occurs when `nextfloat(x0) = x2`.
## The bracket `[a,b]` must be bounded.

function init_state{T <: Float64}(method::Bisection, fs, x::Vector{T}, bracket)
    x0, x2 = sort(x[1:2])
    isinf(x0) && (x0 = nextfloat(x0))
    isinf(x2) && (x2 = prevfloat(x2))
    y0, y2 = fs.f.([x0, x2])

    sign(y0) * sign(y2) > 0 && throw(ArgumentError(bracketing_error))

    state = UnivariateZeroState(x2, x0,
                                y2, y0,
                                isa(bracket, Nullable) ? bracket : Nullable(convert(Vector{T}, sort(bracket))),
                                0, 2,
                                false, false, false, false,
                                "")
    state
end


function update_state{T<:Float64,S}(method::Bisection, fs, o::Roots.UnivariateZeroState{T,S}, options::UnivariateZeroOptions)
    f = fs.f
    x0, x2 = o.xn0, o.xn1
    y0, y2 = o.fxn0, o.fxn1

    x1 = _middle(x0, x2)

    y1 = f(x1)
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

## convergence is much differen there, as we bound between x0, nextfloat(x0) is not measured by eps(), but eps(x0)
function assess_convergence(method::Bisection, fs, state, options)
    x0, x2 = state.xn0, state.xn1
    if x0 > x2
        x0, x2 = x2, x0
    end

    x1 = _middle(x0, x2)

    x0 < x1 && x1 < x2 && return false

    state.message = ""
    state.x_converged=true
    true
end


##################################################

"""

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
                   verbose::Bool=false
             )
    if a > b
        a,b = b,a
    end

    if a >= b || sign(f(a))*sign(f(b)) >= 0
        error("on input a < b and f(a)f(b) < 0 must both hold")
    end
    if xtol < 0.0
        error("tolerance must be >= 0.0")
    end

    c = secant(f, a, b)
    a42a(f, a, b, c,
         xtol=xtol, maxeval=maxeval, verbose=verbose)
end

"""

Split Alefeld, F. A. Potra, and Y. Shi algorithm 4.2 into a function
where `c` is passed in.

Solve f(x) = 0 over bracketing interval [a,b] starting at c, with a < c < b

"""
function a42a(f, a, b, c=(0.5)*(a+b);
                   xtol=zero(float(a)),
                   maxeval::Int=15,
                   verbose::Bool=false
              )


    try
        # re-bracket and check termination
        a, b, d = bracket(f, a, b, c, xtol)
        ee = d
        for n = 2:maxeval
            # use either a cubic (if possible) or quadratic interpolation
            if n > 2 && distinct(f, a, b, d, ee)
                c = ipzero(f, a, b, d, ee)
            else
                c = newton_quadratic(f, a, b, d, 2)
            end
            # re-bracket and check termination
            ab, bb, db = bracket(f, a, b, c, xtol)
            eb = d
            # use another cubic (if possible) or quadratic interpolation
            if distinct(f, ab, bb, db, eb)
                cb = ipzero(f, ab, bb, db, eb)
            else
                cb = newton_quadratic(f, ab, bb, db, 3)
            end
            # re-bracket and check termination
            ab, bb, db = bracket(f, ab, bb, cb, xtol)
            # double length secant step; if we fail, use bisection
            u = abs(f(ab)) < abs(f(bb)) ? ab : bb
            cb = u - 2*f(u)/(f(bb) - f(ab))*(bb - ab)
            ch = abs(cb - u) > (bb - ab)/2 ? ab + (bb - ab)/2 : cb
            # re-bracket and check termination
            ah, bh, dh = bracket(f, ab, bb, ch, xtol)
            # if not converging fast enough bracket again on a bisection
            if bh - ah < 0.5*(b - a)
                a = ah
                b = bh
                d = dh
                ee = db
            else
                ee = dh
                a, b, d = bracket(f, ah, bh, ah + (bh - ah)/2,
                                  xtol)
            end

            verbose && println(@sprintf("a=%18.15f, n=%s", float(a), n))

            if nextfloat(ch) * prevfloat(ch) <= 0
                throw(StateConverged(ch))
            end
            if nextfloat(a) >= b
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
function tole{S,R}(a::S, b::R, fa, fb, tol)
    T = promote_type(S,R)
    u = abs(fa) < abs(fb) ? abs(a) : abs(b)
    2u*eps(T) + tol
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
function bracket(f, a, b, c, tol)

    fa = f(a)
    fb = f(b)

    if !(a <= c <= b)
        error("c must be in (a,b)")
    end
    delta = 0.7*tole(a, b, fa, fb, tol)
    if b - a <= 4delta
        c = (a + b)/2
    elseif c <= a + 2delta
        c = a + 2delta
    elseif c >= b - 2delta
        c = b - 2delta
    end
    fc = f(c)
    if fc == 0
        throw(StateConverged(c))
    elseif sign(fa)*sign(fc) < 0
        aa = a
        bb = c
        db = b
    else
        aa = c
        bb = b
        db = a
    end
    faa = f(aa)
    fbb = f(bb)
    if bb - aa < 2*tole(aa, bb, faa, fbb, tol)
        x0 = abs(faa) < abs(fbb) ? aa : bb
        throw(StateConverged(x0))
    end
    aa, bb, db
end


# take a secant step, if the resulting guess is very close to a or b, then
# use bisection instead
function secant{T}(f, a::T, b)
    c = a - f(a)/(f(b) - f(a))*(b - a)
    tol = eps(T)*5
    if isnan(c) || c <= a + abs(a)*tol || c >= b - abs(b)*tol
        return a + (b - a)/2
    end
    return c
end


# approximate zero of f using quadratic interpolation
# if the new guess is outside [a, b] we use a secant step instead
# based on algorithm on page 330 of [1]
function newton_quadratic(f, a, b, d, k::Int)
    fa = f(a)
    fb = f(b)
    fd = f(d)
    B = (fb - fa)/(b - a)
    A = ((fd - fb)/(d - b) - B)/(d - a)
    if A == 0
        return secant(f, a, b)
    end
    r = A*fa > 0 ? a : b
    for i = 1:k
        r -= (fa + (B + A*(r - b))*(r - a))/(B + A*(2*r - a - b))
    end
    if isnan(r) || (r <= a || r >= b)
        r = secant(f, a, b)
    end
    return r
end


# approximate zero of f using inverse cubic interpolation
# if the new guess is outside [a, b] we use a quadratic step instead
# based on algorithm on page 333 of [1]
function ipzero(f, a, b, c, d)
    fa = f(a)
    fb = f(b)
    fc = f(c)
    fd = f(d)
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
        c = newton_quadratic(f, a, b, d, 3)
    end
    return c
end


# floating point comparison function
function almost_equal(x, y)
    const min_diff = realmin(Float64)*32
    abs(x - y) < min_diff
end


# check that all interpolation values are distinct
function distinct(f, a, b, d, e)
    f1 = f(a)
    f2 = f(b)
    f3 = f(d)
    f4 = f(e)
    !(almost_equal(f1, f2) || almost_equal(f1, f3) || almost_equal(f1, f4) ||
      almost_equal(f2, f3) || almost_equal(f2, f4) || almost_equal(f3, f4))
end


"""
Searches for zeros  of `f` in an interval [a, b].

Basic algorithm used:

* split interval [a,b] into `no_pts` subintervals.
* For each bracketing interval finds a bracketed zero.
* For other subintervals does a quick search with a derivative-free method.

If there are many zeros relative to the number of points, the process
is repeated with more points, in hopes of finding more zeros for
oscillating functions.

"""
function find_zeros(f, a::Real, b::Real, args...;
                    no_pts::Int=100,
                    ftol::Real=10*eps(), reltol::Real=10*eps(),
                    kwargs...)

    a, b = a < b ? (a,b) : (b,a)
    rts = eltype(promote(float(a),b))[]
    xs = vcat(a, a + (b-a)*sort(rand(no_pts)), b)


    ## Look in [ai, bi)
    for i in 1:(no_pts+1)
        ai,bi=xs[i:i+1]
        if isapprox(f(ai), 0.0, rtol=reltol, atol=ftol)
            push!(rts, ai)
        elseif sign(f(ai)) * sign(f(bi)) < 0
            push!(rts, find_zero(f, [ai, bi], Bisection()))
        else
            try
                x = find_zero(f, ai + (0.5)* (bi-ai), Order8(); maxevals=10, reltol=ftol, xreltol=reltol)
                if ai < x < bi
                    push!(rts, x)
                end
            catch e
            end
        end
    end
    ## finally, b?
    isapprox(f(b), 0.0, rtol=reltol, atol=ftol) && push!(rts, b)

    ## redo if it appears function oscillates alot in this interval...
    if length(rts) > (1/4) * no_pts
        return(find_zeros(f, a, b, args...; no_pts = 10*no_pts, ftol=ftol, reltol=reltol, kwargs...))
    else
        return(sort(rts))
    end
end
