"""
    Order16()
    Thukral16()

Implements the order 16 algorithm from
*New Sixteenth-Order Derivative-Free Methods for Solving Nonlinear Equations*
by R. Thukral,
American Journal of Computational and Applied Mathematics
p-ISSN: 2165-8935;    e-ISSN: 2165-8943; 2012;  2(3): 112-118
DOI: [10.5923/j.ajcam.20120203.08](https://doi.org/10.5923/j.ajcam.20120203.08).

Five function calls per step are required. Though rapidly converging,
this method generally isn't faster (fewer function calls/steps) over
other methods when using `Float64` values, but may be useful for
solving over `BigFloat`.  The `Order16` method replaces a Steffensen step with a secant
step when `f(x)` is large.

The error, `eᵢ = xᵢ - α`, is expressed as `eᵢ₊₁ = K⋅eᵢ¹⁶` for an explicit `K`
in equation (50) of the paper.

"""
struct Order16 <: AbstractSecantMethod end
struct Thukral16 <: AbstractSecantMethod end

function update_state(
    M::Order16,
    fs,
    o::UnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    update_state_guarded(M, Secant(), Thukral16(), fs, o, options, l)
end

function update_state(
    M::Thukral16,
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

    fp, issue = _fbracket(xn, wn, fxn, fwn)

    if issue
        log_message(l, "issue with f[xn,wn]")
        o = _set(o, (wn, fwn), (xn, fxn))

        return o, true
    end

    yn::T = xn - fxn / fp
    fyn::S = F(yn)
    incfn(l)

    fp, issue = _fbracket_ratio(yn, xn, wn, fyn, fxn, fwn)
    if issue
        log_message(l, "issue with f[xn,yn]*f[yn,wn]/f[xn,wn]")
        o = _set(o, (yn, fyn), (xn, fxn))

        return o, true
    end

    zn::T = yn - fyn / fp
    fzn::S = F(zn)
    incfn(l)

    fp, issue = _fbracket_diff(xn, yn, zn, fxn, fyn, fzn)
    u2, u3, u4 = fzn / fwn, fyn / fxn, fyn / fwn

    eta = 1 / (1 + 2 * u3 * u4^2) / (1 - u2)
    if issue
        log_message(l, "Approximate derivative failed")
        o = _set(o, (zn, fzn), (xn, fxn))

        return o, true
    end

    an::T = zn - eta * fzn / fp
    fan::S = F(an)
    incfn(l)

    fp, issue = _fbracket_ratio(an, yn, zn, fan, fyn, fzn)
    if issue
        log_message(l, "Approximate derivative failed")
        o = _set(o, (an, fan), (xn, fxn))

        return o, true
    end

    u1, u5, u6 = fzn / fxn, fan / fxn, fan / fwn

    fp1, issue = _fbracket(xn, yn, fxn, fyn)

    sigma =
        1 + u1 * u2 - u1 * u3 * u4^2 +
        u5 +
        u6 +
        u1^2 * u4 +
        u2^2 * u3 +
        3 * u1 * u4^2 * (u3^2 - u4^2) / (fp1 / oneunit(fp1))

    xn1::T = an - sigma * fan / fp
    fxn1::S = F(xn1)
    incfn(l)

    o = _set(o, (xn1, fxn1), (xn, fxn))

    return o, false
end

##################################################
## some means of guarding against large fx when taking a steffensen step
## for Orders 5, 8, and 16 we restrict the size of the steffensen step
## to make the underlying methods more robust. As we need the types defined
## below, this must be included after Order5, Order8, and Order16 are defined
##
function steff_step(M::Union{Order5,Order8,Order16}, x::S, fx::T) where {S,T}
    xbar, fxbar = real(x / oneunit(x)), fx / oneunit(fx)
    thresh = max(1, abs(xbar)) * sqrt(eps(one(xbar))) #^(1/2) # max(1, sqrt(abs(x/fx))) * 1e-6

    out = abs(fxbar) <= thresh ? fxbar : sign(fx) * thresh
    x + out * oneunit(x)
end
