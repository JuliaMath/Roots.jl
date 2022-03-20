"""
    Roots.Brent()

An implementation of
[Brent's](https://en.wikipedia.org/wiki/Brent%27s_method) (or Brent-Dekker) method.
This method uses a choice of inverse quadratic interpolation or a secant
step, falling back on bisection if necessary.

"""
struct Brent <: AbstractAcceleratedBisection end

struct BrentState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    c::T
    d::T
    fxn1::S
    fxn0::S
    fc::S
    mflag::Bool
end


# # we store mflag as -1, or +1 in state.mflag
function init_state(::Brent, F, x₀, x₁, fx₀, fx₁)
    u, v, fu, fv = x₀, x₁, fx₀, fx₁

    if abs(fu) > abs(fv)
        u, v, fu, fv = v, u, fv, fu
    end

    # check if fu*fv ≥ 0
    (iszero(fu) || iszero(fv)) && return BrentState(u, v, v, v, fu, fv, fv, true)
    assert_bracket(fu, fv)

    BrentState(u, v, v, v, fu, fv, fv, true)
end

default_tolerances(::Brent, ::Type{T}, ::Type{S}) where {T,S} = default_tolerances(Secant(), T, S) # need relaxing

function update_state(
    ::Brent,
    f,
    state::BrentState,
    options::UnivariateZeroOptions{T,S},
    l=NullTracks(),
) where {T,S}
    mflag = state.mflag
    a, b, c, d = state.xn0, state.xn1, state.c, state.d
    fa, fb, fc = state.fxn0, state.fxn1, state.fc

    # next step depends on points; inverse quadratic
    s::T = inverse_quadratic_step(a, b, c, fa, fb, fc)
    (isnan(s) || isinf(s)) && (s = secant_step(a, b, fa, fb))

    # guard step
    u, v = (3a + b) / 4, b
    if u > v
        u, v = v, u
    end

    tol = max(options.xabstol, max(abs(b), abs(c), abs(d)) * options.xreltol)
    if !(u < s < v) ||
       (mflag && abs(s - b) >= abs(b - c) / 2) ||
       (!mflag && abs(s - b) >= abs(b - c) / 2) ||
       (mflag && abs(b - c) <= tol) ||
       (!mflag && abs(c - d) <= tol)
        s = _middle(a, b)

        mflag = true
    else
        mflag = false
    end

    fs = f(s)
    incfn(l)

    if iszero(fs)
        @set! state.xn1 = s
        @set! state.fxn1 = fs
        return (state, true)
    elseif isnan(fs) || isinf(fs)
        return (state, true)
    end

    d = c
    c, fc = b, fb

    if sign(fa) * sign(fs) < 0
        b, fb = s, fs
    else
        a, fa = s, fs
    end

    if abs(fa) < abs(fb)
        a, b, fa, fb = b, a, fb, fa
    end

    @set! state.xn0 = a
    @set! state.xn1 = b
    @set! state.c = c
    @set! state.d = d
    @set! state.fxn0 = fa
    @set! state.fxn1 = fb
    @set! state.fc = fc
    @set! state.mflag = mflag

    return state, false
end
