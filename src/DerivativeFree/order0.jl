"""
    Order0()


The `Order0` method is engineered to be a more robust, though possibly
slower, alternative to the other derivative-free root-finding
methods. The implementation roughly follows the algorithm described in
*Personal Calculator Has Key to Solve Any Equation ``f(x) = 0``*, the
SOLVE button from the
[HP-34C](http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf).
The basic idea is to use a secant step. If along the way a bracket is
found, switch to a bracketing algorithm, using `AlefeldPotraShi`.  If the secant
step fails to decrease the function value, a quadratic step is used up
to ``4`` times.

This is not really ``0``-order: the secant method has order
``1.6...`` [Wikipedia](https://en.wikipedia.org/wiki/Secant_method#Comparison_with_other_root-finding_methods)
and the the bracketing method has order
``1.6180...`` [Wikipedia](http://www.ams.org/journals/mcom/1993-61-204/S0025-5718-1993-1192965-2/S0025-5718-1993-1192965-2.pdf)
so for reasonable starting points and functions, this algorithm should be
superlinear, and relatively robust to non-reasonable starting points.

"""
struct Order0 <: AbstractSecantMethod end

# special case Order0 to be hybrid
function init(
    𝑭𝑿::ZeroProblem,
    M::Order0,
    p′=nothing;
    p=nothing,
    verbose::Bool=false,
    tracks=NullTracks(),
    kwargs...,
)
    p = p′ === nothing ? p : p′
    init(𝑭𝑿, Secant(), AlefeldPotraShi(); p=p, verbose=verbose, tracks=tracks, kwargs...)
end

init(::ZeroProblem, ::Order0, ::AbstractBracketingMethod; kwargs...) =
    throw(ArgumentError("No bracketing method specified with Order0"))
