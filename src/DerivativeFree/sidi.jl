"""
    Sidi(k)

Implements family of methods in "Generalization Of The Secant Method For Nonlinear Equations" by  Avram Sidi in Applied Mathematics E-Notes, 8(2008), 115-123.

These methods generalize the secant method by using an interpolating polynomial based on the last ``k+1`` points to estimate ``f'(xₙ)`` in its use with Newton's method, the secant method being the ``k=1`` case.


## Convergence rates:

* `Sidi(1) = 1.618⋯` (secant method)
* `Sidi(2) = 1.839⋯`
* `Sidi(3) = 1.928⋯`
* `Sidi(4) = 1.966⋯`

Like the secant method, each update step requires one function evaluation.

## Example

```
find_zero(sin, 3, Roots.Sidi(2))
```

"""
struct Sidi{k} <: AbstractSecantMethod end
Sidi(k::Int) = Sidi{k}()

struct SidiState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    fxn1::S
    fxn0::S
    xs::Vector{T}
    fs::Vector{T}
end

function init_state(M::Sidi{k}, F::Callable_Function, x) where {k}
    x₀, x₁ = x₀x₁(x)
    fx₀, xs, fs = _init_sidi(F, (x₀, x₁), k)
    state = SidiState(xs[k], xs[k + 1], fx₀, fs[1], xs, fs)
end

function update_state(
    L::Sidi{k},
    F,
    o::SidiState{T,S},
    options,
    l=NullTracks(),
) where {k,T,S}
    xs, fs = o.xs, o.fs
    fxn1 = o.fxn1
    _update_sidi!(F, xs, fs)

    incfn(l)
    @reset o.xn0 = xs[k]
    @reset o.xn1 = xs[k + 1]
    @reset o.fxn0 = o.fxn1
    @reset o.fxn1 = fs[1]
    @reset o.xs = xs
    @reset o.fs = fs

    return (o, false)
end

# create xs and fs; fs records lower diagonal in newton tableau
# x1  .  .    .  f1234
# x2  .  .  f234
# x3  . f34
# x4 f4
# xs = [x1,x2,x3,x4]; fs = [f4, f34, f234, f1234]
# build up diagonal by diagonal
# where pk' uses xs, fs to be computed
# no good way to pre-compute fs to speed up, so here
# we expect x to have 2 or more elements
# but we compute each f(x)
# This allocates, as it uses a vector to store xs, fs
function _init_sidi(f, x, k)
    x₀ = first(x)
    fx₀ = f(x₀)

    xs = Vector{typeof(x₀)}(undef, k+1)
    fs = Vector{typeof(fx₀)}(undef, k+1)

    n = length(x)
    xs[1:n] .= x

    xs[1] = x[1]
    xs[2] = x[2]

    # diagonal for x2 above (two entries)
    fs[1] = f(xs[2])
    fs[2] = (fx₀ - fs[1]) / (xs[1] - xs[2])

    # build up diagonal by diagonal
    for j in 3:(k + 1)
        if j ≤ n # xⱼ was specified
            xⱼ = xs[j]
        else
            xⱼ₋₁ = xs[j - 1]
            pk′ = evaluate_pk′(view(xs, 1:(j - 1)), view(fs, 1:(j - 1)))
            xⱼ = xs[j] = xⱼ₋₁ - fs[1] / pk′
        end
        Δ = f(xⱼ)
        for i in 2:j
            Δ₀ = fs[i - 1]
            fs[i - 1] = Δ
            Δ = (Δ₀ - Δ) / (xs[j - i + 1] - xs[j])
        end
        fs[j] = Δ
    end

    # return fx₀ for bookkeeping purposes
    fx₀, xs, fs
end

# update step: compute xn, fxn, update the xs,fs tableau
function _update_sidi!(f, xs, fs)
    xₙ₋₁, fxₙ₋₁ = xs[end], fs[1]
    fn′ = evaluate_pk′(xs, fs)
    xn = xₙ₋₁ - fxₙ₋₁ / fn′
    fxn = f(xn)
    update_tableau!(xn, fxn, xs, fs)
end

# formula (10) in paper to evaluate derivative of interpolating polynomial
function evaluate_pk′(xs1, fs1)
    δ = xs1[end] - xs1[end - 1]
    Σ = fs1[2]
    k = length(xs1)
    for i in 3:k
        Σ = Σ + fs1[i] * δ
        δ = δ * (xs1[end] - xs1[end - i + 1])
    end
    Σ
end

# update tableau's lower part
# leaves [xn-k, xn-k+1, xn-k+2, ..., xn]
#        [fn, f(n-1,n), f(n-2, n-1, n), ..., f(n-k,n-k+1, ..., n)]
function update_tableau!(xn, fxn, xs0, fs0)
    k = length(xs0)
    for i in 1:(k - 1)
        xs0[i] = xs0[i + 1]
    end
    xs0[end] = xn
    Δ = fxn
    for i in 2:k
        Δ₀ = fs0[i - 1]
        fs0[i - 1] = Δ
        Δ = (Δ₀ - Δ) / (xs0[end - i + 1] - xn)
    end
    fs0[end] = Δ
    xs0, fs0
end
