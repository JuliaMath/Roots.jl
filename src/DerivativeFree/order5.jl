"""
    Order5()
    KumarSinghAkanksha()

Implements an order 5 algorithm from *A New Fifth Order Derivative
Free Newton-Type Method for Solving Nonlinear Equations* by Manoj
Kumar, Akhilesh Kumar Singh, and Akanksha, Appl. Math. Inf. Sci. 9,
No. 3, 1507-1513 (2015), DOI: [10.12785/amis/090346](https://doi.org/10.12785/amis/090346). Four function
calls per step are needed.  The `Order5` method replaces a Steffensen step with a secant
step when `f(x)` is large.

The error, `eᵢ = xᵢ - α`, satisfies
`eᵢ₊₁ = K₁ ⋅ K₅ ⋅ M ⋅ eᵢ⁵ + O(eᵢ⁶)`

"""
struct Order5 <: AbstractSecantMethod end
struct KumarSinghAkanksha <: AbstractSecantMethod end

function update_state(M::Order5, fs, o::UnivariateZeroState, options, l=NullTracks())
    update_state_guarded(M, Secant(), KumarSinghAkanksha(), fs, o, options, l)
end

function update_state(
    M::KumarSinghAkanksha,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    xn = o.xn1
    fxn = o.fxn1

    wn::T = steff_step(M, o.xn1, o.fxn1)

    fwn::S = F(wn)
    incfn(l)

    fp, issue = _fbracket(o.xn1, wn, o.fxn1, fwn)
    if issue
        log_message(l, "Issue with divided difference f[xn, wn]. ")
        o = _set(o, (wn, fwn), (xn, fxn))
        return o, true
    end

    yn::T = o.xn1 - o.fxn1 / fp
    fyn::S = F(yn)
    incfn(l)

    zn::T = xn - (fxn + fyn) / fp
    fzn::S = F(zn)
    incfn(l)

    fp, issue = _fbracket_ratio(yn, o.xn1, wn, fyn, o.fxn1, fwn)
    if issue
        log_message(l, "Issue with f[xn,yn] * f[yn,wn] / f[xn, wn].")
        o = _set(o, (yn, fyn), (xn, fxn))

        return o, true
    end

    x₁::T = zn - fzn / fp
    f₁ = F(x₁)
    incfn(l)

    o = _set(o, (x₁, f₁), (xn, fxn))

    return o, false

    #    nothing
end
