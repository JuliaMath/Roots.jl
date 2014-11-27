## Bisection_method for floats and int

# type to throw on succesful convergence
type StateConverged
    x0::Real
end


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
 
  xint = reinterpret(Uint64, abs(x))
  yint = reinterpret(Uint64, abs(y))
  unsigned = reinterpret(Float64, (xint + yint) >> 1)
 
  return negate ? -unsigned : unsigned
end
 

function find_zero(f::Function, a::Float64, b::Float64)

    x0, y0 = a, f(a)
    x2, y2 = b, f(b)

    y0 == 0 && return a
    y2 == 0.0 && return b
    if sign(y0) * sign(y2) > 0
        error("[a,b] is not a bracket") 
    end

    x1 = _middle(x0, x2)
    y1 = f(x1)

    
    while x0 < x1 && x1 < x2
        if sign(y0) == sign(y1)
            x0, x2 = x1, x2
            y0, y2 = y1, y2
        else
            x0, x2 = x0, x1
            y0, y2 = y0, y1
        end
        
        x1 = _middle(x0, x2)
        y1 = f(x1)
        sign(y1) == 0 && return x1
    end
    
    return abs(y0) < abs(y2) ? x0 : x2
end


# fzero finds the root of a continuous function within a provided
# interval [a, b], without requiring derivatives. It is based on algorithm 4.2
# described in: 1. G. E. Alefeld, F. A. Potra, and Y. Shi, "Algorithm 748:
# enclosing zeros of continuous functions," ACM Trans. Math. Softw. 21,
# 327–344 (1995).
#
# input:
#     f: function to find the root of
#     a, b: the initial bracket, with: a < b, f(a)*f(b) < 0
#     tolerance: acceptable error (it's safe to set zero for machine precision)
#     max_iter:  maximum number of iterations
#
# output:
#     an estimate of the zero of f
function find_zero(f::Function, a::Union(BigFloat, BigInt), b::Union(BigFloat, BigInt);
                   tol=0.0, 
                   max_iter=100,
                   verbose::Bool=false
                   )
    if a >= b || sign(f(a))*sign(f(b)) >= 0
        error("on input a < b and f(a)f(b) < 0 must both hold")
    end
    if tol < 0.0
        error("tolerance must be >= 0.0")
    end
    try
        # start with a secant approximation
        c = secant(f, a, b)
        # re-bracket and check termination
        a, b, d = bracket(f, a, b, c, tol)
        for n = 2:max_iter
            # use either a cubic (if possible) or quadratic interpolation
            if n > 2 && distinct(f, a, b, d, e)
                c = ipzero(f, a, b, d, e)
            else
                c = newton_quadratic(f, a, b, d, 2)
            end
            # re-bracket and check termination
            ab, bb, db = bracket(f, a, b, c, tol)
            eb = d
            # use another cubic (if possible) or quadratic interpolation
            if distinct(f, ab, bb, db, eb)
                cb = ipzero(f, ab, bb, db, eb)
            else
                cb = newton_quadratic(f, ab, bb, db, 3)
            end
            # re-bracket and check termination
            ab, bb, db = bracket(f, ab, bb, cb, tol)
            # double length secant step; if we fail, use bisection
            u = abs(f(ab)) < abs(f(bb)) ? ab : bb
            cb = u - 2*f(u)/(f(bb) - f(ab))*(bb - ab)
            ch = abs(cb - u) > (bb - ab)/2 ? ab + (bb - ab)/2 : cb
            # re-bracket and check termination
            ah, bh, dh = bracket(f, ab, bb, ch, tol)
            # if not converging fast enough bracket again on a bisection
            if bh - ah < 0.5*(b - a)
                a = ah
                b = bh
                d = dh
                e = db
            else
                e = dh
                a, b, d = bracket(f, ah, bh, ah + (bh - ah)/2,
                                  tol)
            end

            verbose && println("a=$a, n=$n")
        end
    catch ex
        if isa(ex, StateConverged)
            return ex.x0
        else
            throw(ex)
        end
    end
    error("root not found within max_iter iterations")
end



# type to throw on failure
type ConvergenceFailed end

# calculate a scaled tolerance
# based on algorithm on page 340 of [1]
function tole(a, b, fa, fb, tol)
    u = abs(fa) < abs(fb) ? abs(a) : abs(b)
    2u*eps(1.0) + tol
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
function bracket(f::Function, a, b, c, tol)
    if !(a <= c <= b)
        error("c must be in (a,b)")
    end
    fa = f(a)
    fb = f(b)
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
function secant(f::Function, a, b)
    c = a - f(a)/(f(b) - f(a))*(b - a)
    tol = eps(1.0)*5
    if c <= a + abs(a)*tol || c >= b - abs(b)*tol
        return a + (b - a)/2
    end
    return c
end


# approximate zero of f using quadratic interpolation
# if the new guess is outside [a, b] we use a secant step instead
# based on algorithm on page 330 of [1]
function newton_quadratic(f::Function, a, b, d, k::Int)
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
    if r <= a || r >= b
        r = secant(f, a, b)
    end
    return r
end


# approximate zero of f using inverse cubic interpolation
# if the new guess is outside [a, b] we use a quadratic step instead 
# based on algorithm on page 333 of [1]
function ipzero(f::Function, a, b, c, d)
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
function distinct(f::Function, a, b, d, e)
    f1 = f(a)
    f2 = f(b)
    f3 = f(d)
    f4 = f(e)
    !(almost_equal(f1, f2) || almost_equal(f1, f3) || almost_equal(f1, f4) ||
      almost_equal(f2, f3) || almost_equal(f2, f4) || almost_equal(f3, f4))
end


## split interval [a,b] into no_pts intervals, apply find_zero to each, accumulae
function find_zeros(f::Function, a::Real, b::Real, args...;no_pts::Int=200, kwargs...)
    a, b = a < b ? (a,b) : (b,a)

    rts = eltype(promote(float(a),b))[]
    xs = linspace(a, b, no_pts)    

    for i in 1:(length(xs)-1)
        if f(xs[i]) * f(xs[i+1]) < 0
            push!(rts, fzero(f, xs[i:(i+1)]))
        end
        if f(xs[i]) == 0.0
            push!(rts, xs[i])
        end
    end
    if f(b) == 0
        push!(rts, b)
    end

    rts
end
