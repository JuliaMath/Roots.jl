## --------------------------------------------------
"""
    find_bracket(f, x0, method=A42(); kwargs...)

For bracketing methods returns an approximate root, the last bracketing interval used, and a flag indicating if an exact zero was found as a named tuple.

With the default tolerances, one  of these should be the case: `exact`  is `true` (indicating termination  of the algorithm due to an exact zero  being identified) or the length of `bracket` is less or equal than `2eps(maximum(abs.(bracket)))`. In the `BisectionExact` case, the 2 could  be replaced by 1, as the bracket, `(a,b)` will satisfy  `nextfloat(a) == b `;  the Alefeld,  Potra, and Shi algorithms don't quite have  that promise.

"""
function find_bracket(
    fs,
    x0,
    method::M=A42();
    kwargs...,
) where {M<:Union{AbstractAlefeldPotraShi,BisectionExact}}

    Base.depwarn("This interface is deprecated", :find_bracket)

    x = adjust_bracket(x0)
    F = Callable_Function(method, fs) #callable_function(fs)
    state = init_state(method, F, x)
    options = init_options(method, state; kwargs...)
    l = NullTracks()

    # check if tolerances are exactly 0
    iszero_tol =
        iszero(options.xabstol) &&
        iszero(options.xreltol) &&
        iszero(options.abstol) &&
        iszero(options.reltol)

    val, stopped = :not_converged, false
    while !stopped
        val, stopped = assess_convergence(method, state, options)
        stopped && break
        state, stopped = update_state(method, F, state, options, l)
    end

    a, b = state.xn0, state.xn1
    fa, fb = state.fxn0, state.fxn1
    xstar, fxstar = abs(fa) < abs(fb) ? (a, fa) : (b, fb)
    (xstar=xstar, bracket=(a, b), exact=iszero(fxstar))
end
