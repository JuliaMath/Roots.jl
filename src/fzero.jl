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
function fzero(f::Function, a, b;
               tolerance=0.0, max_iter=100)
    if a >= b || sign(f(a))*sign(f(b)) >= 0
        error("on input a < b and f(a)f(b) < 0 must both hold")
    end
    if tolerance < 0.0
        error("tolerance must be >= 0.0")
    end
    try
        # start with a secant approximation
        c = secant(f, a, b)
        # re-bracket and check termination
        a, b, d = bracket(f, a, b, c, tolerance)
        for n = 2:max_iter
            # use either a cubic (if possible) or quadratic interpolation
            if n > 2 && distinct(f, a, b, d, e)
                c = ipzero(f, a, b, d, e)
            else
                c = newton_quadratic(f, a, b, d, 2)
            end
            # re-bracket and check termination
            ab, bb, db = bracket(f, a, b, c, tolerance)
            eb = d
            # use another cubic (if possible) or quadratic interpolation
            if distinct(f, ab, bb, db, eb)
                cb = ipzero(f, ab, bb, db, eb)
            else
                cb = newton_quadratic(f, ab, bb, db, 3)
            end
            # re-bracket and check termination
            ab, bb, db = bracket(f, ab, bb, cb, tolerance)
            # double length secant step; if we fail, use bisection
            u = abs(f(ab)) < abs(f(bb)) ? ab : bb
            cb = u - 2*f(u)/(f(bb) - f(ab))*(bb - ab)
            ch = abs(cb - u) > (bb - ab)/2 ? ab + (bb - ab)/2 : cb
            # re-bracket and check termination
            ah, bh, dh = bracket(f, ab, bb, ch, tolerance)
            # if not converging fast enough bracket again on a bisection
            if bh - ah < 0.5*(b - a)
                a = ah
                b = bh
                d = dh
                e = db
            else
                e = dh
                a, b, d = bracket(f, ah, bh, ah + (bh - ah)/2,
                                  tolerance)
            end
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


# type to throw on succesful convergence
type StateConverged
    x0::Float64
end


# calculate a scaled tolerance
# based on algorithm on page 340 of [1]
function tole(a, b, fa, fb, tolerance)
    u = abs(fa) < abs(fb) ? abs(a) : abs(b)
    2u*eps(1.0) + tolerance
end


# bracket the root
# inputs:
#     - f: the function
#     - a, b: the current bracket with a < b, f(a)f(b) < 0
#     - c within (a,b): current best guess of the root
#     - tolerance: desired accuracy
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
function bracket(f::Function, a, b, c, tolerance)
    if !(a <= c <= b)
        error("c must be in (a,b)")
    end
    fa = f(a)
    fb = f(b)
    delta = 0.7*tole(a, b, fa, fb, tolerance)
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
    if bb - aa < 2*tole(aa, bb, faa, fbb, tolerance)
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

