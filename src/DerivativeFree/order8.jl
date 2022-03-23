"""
    Order8()
    Thukral8()

Implements an eighth-order algorithm from *New Eighth-Order
Derivative-Free Methods for Solving Nonlinear Equations* by Rajinder
Thukral, International Journal of Mathematics and Mathematical
Sciences Volume 2012 (2012), Article ID 493456, 12 pages DOI:
[10.1155/2012/493456](https://doi.org/10.1155/2012/493456). Four
function calls per step are required.  The `Order8` method replaces a
Steffensen step with a secant step when `f(x)` is large.

The error, `eᵢ = xᵢ - α`, is expressed as `eᵢ₊₁ = K ⋅ eᵢ⁸` in
(2.25) of the paper for an explicit `K`.

"""
struct Order8 <: AbstractSecantMethod end
struct Thukral8 <: AbstractSecantMethod end

## cf also: https://doi.org/10.1515/tmj-2017-0049
function update_state(
    M::Order8,
    fs,
    o::UnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    update_state_guarded(M, Secant(), Thukral8(), fs, o, options, l)
end

function update_state(
    M::Thukral8,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    xn = o.xn1
    fxn = o.fxn1

    wn::T = steff_step(M, xn, fxn)
    fwn::S = F(wn)
    incfn(l)

    if isissue(fwn)
        log_message(l, "issue with Steffensen step fwn")
        o = _set(o, (wn, fwn), (xn, fxn))

        return o, true
    end

    fp, issue = _fbracket(xn, wn, fxn, fwn)

    if issue
        log_message(l, "issue with divided difference f[xn, wn]. ")
        return o, true
    end

    yn::T = xn - fxn / fp
    fyn::S = F(yn)
    incfn(l)

    fp, issue = _fbracket(yn, xn, fyn, fxn)
    if issue #fp
        log_message(l, "issue with divided difference f[xn, yn]. ")
        o = _set(o, (yn, fyn), (xn, fxn))

        return o, true
    end

    phi = (1 + fyn / fwn)           # pick one of options
    zn::T = yn - phi * fyn / fp
    fzn::S = F(zn)
    incfn(l)

    fp, issue = _fbracket_diff(xn, yn, zn, fxn, fyn, fzn)
    if issue
        log_message(l, "issue with divided difference  f[y,z] - f[x,y] + f[x,z]. ")
        o = _set(o, (zn, fzn), (xn, fxn))

        return o, true
    end

    w = 1 / (1 - fzn / fwn)

    xi = (1 - 2fyn * fyn * fyn / (fwn * fwn * fxn))

    xn1::T = zn - w * xi * fzn / fp
    fxn1::S = F(xn1)
    incfn(l)

    o = _set(o, (xn1, fxn1), (xn, fxn))

    return o, false
end
