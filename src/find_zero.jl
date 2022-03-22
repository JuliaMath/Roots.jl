"""

    find_zero(f, x0, M, [N::AbstractBracketing]; kwargs...)

Interface to one of several methods for finding zeros of a univariate function, e.g. solving ``f(x)=0``.

# Initial starting value

For most methods, `x0` is a scalar value indicating the initial value
in the iterative procedure. (Secant methods can have a tuple specify
their initial values.) Values must be a subtype of `Number` and have
methods for `float`, `real`, and `oneunit` defined.

For bracketing intervals, `x0` is specified using a tuple, a vector, or any iterable with `extrema` defined. A bracketing interval, ``[a,b]``, is one where f(a) and f(b) have different signs.

# Return value

If the algorithm suceeds, the approximate root identified is returned. A `ConvergenceFailed` error is thrown if the algorithm fails. The alternate form `solve(ZeroProblem(f,x0), M)` returns `NaN` in case of failure.

# Specifying a method

A method is specified to indicate which algorithm to employ:

* There are methods for bisection where a bracket is specified: [`Bisection`](@ref), [`Roots.A42`](@ref), [`Roots.AlefeldPotraShi`](@ref), [`Roots.Brent`](@ref), and [`FalsePosition`](@ref)

* There are several derivative-free methods: cf. [`Order0`](@ref), [`Order1`](@ref) (also [`Roots.Secant`](@ref)), [`Order2`](@ref) (also [`Roots.Steffensen`](@ref)), [`Order5`](@ref), [`Order8`](@ref), and [`Order16`](@ref), where the number indicates the order of the convergence. Methods [`Roots.Order1B`](@ref) and [`Roots.Order2B`](@ref) are useful when the desired zero has a multiplicity.

* There are some classical methods where derivatives are required: [`Roots.Newton`](@ref), [`Roots.Halley`](@ref), [`Roots.Schroder`](@ref).

* The family [`Roots.LithBoonkkampIJzerman{S,D}`](@ref) for different `S` and `D` uses a linear multistep method root finder. The `(2,0)` method is the secant method, `(1,1)` is Newton's method.

For more detail, see the help page for each method (e.g., `?Order1`). Non-exported methods must be qualified with module name, as in `?Roots.Schroder`.

If no method is specified, the default method depends on `x0`:

* If `x0` is a scalar, the default is the slower, but more robust `Order0` method.

* If `x0` is a tuple, vector, or iterable with `extrema` defined indicating a *bracketing* interval, the `Bisection` method is used. (The exact algorithm depends on the number type and the tolerances.)

# Specifying the function

The function(s) are passed as the first argument.

For the few methods that use one or more derivatives (`Newton`, `Halley`,
`Schroder`, `LithBoonkkampIJzerman(S,D)`, and  `Order5Derivative`) a
tuple of functions is used. For the classical algorithms, a function returning `(f(x), f(x)/f'(x), [f'(x)/f''(x)])` may be used.

# Optional arguments (tolerances, limit evaluations, tracing)

* `xatol` - absolute tolerance for `x` values. Passed to `isapprox(x_n, x_{n-1})`.
* `xrtol` - relative tolerance for `x` values. Passed to `isapprox(x_n, x_{n-1})`.
* `atol`  - absolute tolerance for `f(x)` values.
* `rtol`  - relative tolerance for `f(x)` values.
* `maxevals`   - limit on maximum number of iterations.
* `strict` - if `false` (the default), when the algorithm stops, possible zeros are checked with a relaxed tolerance.
* `verbose` - if `true` a trace of the algorithm will be shown on successful completion. See the internal [`Tracks`](@ref) object to save this trace.

See the help string for `Roots.assess_convergence` for details on
convergence. See the help page for `Roots.default_tolerances(method)`
for details on the default tolerances.

In general, with floating point numbers, convergence must be
understood as not an absolute statement. Even if mathematically `α` is
an answer and `xstar` the floating point realization, it may be that
`f(xstar) - f(α)  ≈ xstar ⋅  f'(α) ⋅ eps(α)`, so tolerances must be
appreciated, and at times specified.

For the `Bisection` methods, convergence is guaranteed, so the tolerances are set to be ``0`` by default.

If a bracketing method is passed in after the method specification,
then whenever a bracket is identified during the algorithm, the method
will switch to the bracketing method to identify the zero. (Bracketing
methods are mathematically guaranteed to converge, non-bracketing methods may or may not converge.)
This is what `Order0` does by default, with an initial secant method switching
to the `AlefeldPotraShi` method should a bracket be encountered.

Note: The order of the method is hinted at in the naming scheme. A
scheme is order `r` if, with `eᵢ = xᵢ - α`, `eᵢ₊₁ = C⋅eᵢʳ`. If the
error `eᵢ` is small enough, then essentially the error will gain `r`
times as many leading zeros each step. However, if the error is not
small, this will not be the case. Without good initial guesses, a high
order method may still converge slowly, if at all. The `OrderN`
methods have some heuristics employed to ensure a wider range for
convergence at the cost of not faithfully implementing the method,
though those are available through unexported methods.

# Examples:

Default methods.

```jldoctest find_zero
julia> using Roots

julia> find_zero(sin, 3)  # use Order0()
3.141592653589793

julia> find_zero(sin, (3,4)) # use Bisection()
3.141592653589793
```

Specifying a method,

```jldoctest find_zero
julia> find_zero(sin, (3,4), Order1())            # can specify two starting points for secant method
3.141592653589793

julia> find_zero(sin, 3.0, Order2())              # Use Steffensen method
3.1415926535897936

julia> find_zero(sin, big(3.0), Order16())        # rapid convergence
3.141592653589793238462643383279502884197169399375105820974944592307816406286198

julia> find_zero(sin, (3, 4), Roots.A42())      # fewer function calls than Bisection(), in this case
3.141592653589793

julia> find_zero(sin, (3, 4), FalsePosition(8))   # 1 of 12 possible algorithms for false position
3.141592653589793

julia> find_zero((sin,cos), 3.0, Roots.Newton())  # use Newton's method
3.141592653589793

julia> find_zero((sin, cos, x->-sin(x)), 3.0, Roots.Halley())  # use Halley's method
3.141592653589793
```

Changing tolerances.

```jldoctest find_zero
julia> fn = x -> (2x*cos(x) + x^2 - 3)^10/(x^2 + 1);

julia> x0, xstar = 3.0,  2.9947567209477;

julia> fn(find_zero(fn, x0, Order2())) <= 1e-14  # f(xₙ) ≈ 0, but Δxₙ can be largish
true

julia> find_zero(fn, x0, Order2(), atol=0.0, rtol=0.0) # error: x_n ≉ x_{n-1}; just f(x_n) ≈ 0
ERROR: Roots.ConvergenceFailed("Algorithm failed to converge")
[...]

julia> fn = x -> (sin(x)*cos(x) - x^3 + 1)^9;

julia> x0, xstar = 1.0,  1.112243913023029;

julia> find_zero(fn, x0, Order2()) ≈ xstar
true

julia> find_zero(fn, x0, Order2(), maxevals=3)    # need more steps to converge
ERROR: Roots.ConvergenceFailed("Algorithm failed to converge")
[...]
```

# Tracing

Passing `verbose=true` will show details on the steps of the algorithm.
The `tracks` argument allows
the passing of storage to record the values of `x` and `f(x)` used in
the algorithm.

!!! note
    See [`solve!`](@ref) and [`ZeroProblem`](@ref) for an alternate interface.
"""
function find_zero(
    f,
    x0,
    M::AbstractUnivariateZeroMethod;
    p=nothing,
    verbose=false,
    tracks::AbstractTracks=NullTracks(),
    kwargs...,
)
    xstar = solve(ZeroProblem(f, x0), M; p=p, verbose=verbose, tracks=tracks, kwargs...)

    isnan(xstar) && throw(ConvergenceFailed("Algorithm failed to converge"))

    xstar
end

# defaults when method is not specified
# if a number, use Order0
# O/w use a bracketing method of an assumed iterable
find_zero(f, x0::Number; kwargs...) = find_zero(f, x0, Order0(); kwargs...)
find_zero(f, x0; kwargs...) = find_zero(f, x0, Bisection(); kwargs...)

"""
    find_zero(M, F, state, [options], [l])

Find zero using method `M`, function(s) `F`, and initial state
`state`. Returns an approximate zero or `NaN`. Useful when some part
of the processing pipeline is to be adjusted.

* `M::AbstractUnivariateZeroMethod` a method, such as `Secant()`
* `F`: A callable object (or tuple of callable objects for certain methods)
* `state`: An initial state, as created by `init_state` (or `_init_state`).
* `options::UnivariateZeroOptions`: specification of tolerances
* `l::AbstractTracks`: used to record steps in algorithm, when requested.
```

!!! note
    To be deprecated in favor of `solve!(init(...))`.
"""
function find_zero(
    M::AbstractUnivariateZeroMethod,
    F,
    state::AbstractUnivariateZeroState,
    options::UnivariateZeroOptions=init_options(M, state),
    l::AbstractTracks=NullTracks(),
)
    solve!(init(M, F, state, options, l))
end

## ---------------

## Create an Iterator interface
# returns NaN, not an error, if there are issues

"""
    ZeroProblem{F,X}

A container for a function and initial guess to be used with `solve`.
"""
struct ZeroProblem{F,X}
    F::F
    x₀::X
end
Base.broadcastable(p::ZeroProblem) = Ref(p)

## The actual iterating object
struct ZeroProblemIterator{M,N,F,S,O,L}
    M::M
    N::N
    F::F
    state::S
    options::O
    logger::L
end

Base.show(io::IO, Z::ZeroProblemIterator) =
    print(io, "A problem object to pass to `solve!`")

## Initialize a Zero Problem Iterator
## init(Z,p)
## init(Z,M,p)
## init(M,F,state, [options], [logger])
## want p to be positional, not named⁺
function init(
    𝑭𝑿::ZeroProblem,
    M::AbstractUnivariateZeroMethod,
    p′=nothing;
    p=nothing,
    verbose::Bool=false,
    tracks=NullTracks(),
    kwargs...,
)
    F = Callable_Function(M, 𝑭𝑿.F, p === nothing ? p′ : p)  #⁺
    state = init_state(M, F, 𝑭𝑿.x₀)
    options = init_options(M, state; kwargs...)
    l = Tracks(verbose, tracks, state)
    incfn(l, initial_fncalls(M))
    ZeroProblemIterator(M, nothing, F, state, options, l)
end

function init(𝑭𝑿::ZeroProblem, p′=nothing; kwargs...)
    M = length(𝑭𝑿.x₀) == 1 ? Secant() : Bisection()
    init(𝑭𝑿, M, p′; kwargs...)
end

function init(
    M::AbstractUnivariateZeroMethod,
    F,
    state::AbstractUnivariateZeroState,
    options::UnivariateZeroOptions=init_options(M, state),
    l::AbstractTracks=NullTracks(),
)
    ZeroProblemIterator(M, Nothing, Callable_Function(M, F), state, options, l)
end

# Optional iteration interface to handle looping
# * returns xₙ or (aₙ, bₙ) depending
# * throws error on non-convergence
function Base.iterate(P::ZeroProblemIterator, st=nothing)
    ## st = (val, (state, ctr, flag, stopped))

    M, F, options, l = P.M, P.F, P.options, P.logger

    if st === nothing
        state = P.state
        ctr, flag, stopped = 1, :not_converged, false
        log_method(l, M)
        log_step(l, M, state; init=true)
    else
        state, ctr, flag, stopped = st
        ctr += 1
    end

    stopped && return nothing
    ctr > options.maxevals && return nothing

    state, stopped = update_state(M, F, state, options, l)
    log_step(l, M, state)

    flag, stopped = assess_convergence(M, state, options)

    if stopped
        α = decide_convergence(M, F, state, options, flag)
        log_last(l, α)
        isnan(α) && throw(ConvergenceFailed("Algorithm did not converge."))
    end

    return (last(state,M), (state, ctr, flag, stopped))

end



"""
    solve(fx::ZeroProblem, [p=nothing]; kwargs...)
    solve(fx::ZeroProblem, M, [p=nothing]; kwargs...)
    init(fx::ZeroProblem, [M], [p=nothing];
         verbose=false, tracks=NullTracks(), kwargs...)
    solve!(P::ZeroProblemIterator)

Solve for the zero of a function specified through a  `ZeroProblem` or `ZeroProblemIterator`

The methods involved with this interface are:

* `ZeroProblem`: used to specify a problem with a function (or functions) and an initial guess
* `solve`: to solve for a zero in a `ZeroProblem`

The latter calls the following, which can be useful independently:

* `init`: to initialize an iterator with a method for solution, any adjustments to the default tolerances, and a specification to log the steps or not.
* `solve!` to iterate to convergence.

Returns `NaN`, not an error, when the problem can not be solved. Tested for zero allocations.

## Examples:

```jldoctest find_zero
julia> using Roots

julia> fx = ZeroProblem(sin, 3)
ZeroProblem{typeof(sin), Int64}(sin, 3)

julia> solve(fx)
3.141592653589793
```

Or, if the iterable is required

```jldoctest find_zero
julia> problem = init(fx);


julia> solve!(problem)
3.141592653589793
```

The default method is `Order1()`, when  `x0` is a number, or `Bisection()` when `x0` is an iterable with 2 or more values.


A second position argument for `solve` or `init` is used to specify a different method; keyword arguments can be used to adjust the default tolerances.


```jldoctest find_zero
julia> solve(fx, Order5(), atol=1/100)
3.1425464815525403
```

The above is equivalent to:

```jldoctest find_zero
julia> problem = init(fx, Order5(), atol=1/100);


julia> solve!(problem)
3.1425464815525403
```

The  argument `p` may be used if the function(s) to be solved depend on a parameter in their second positional argument (e.g., `f(x,p)`). For example

```jldoctest find_zero
julia> f(x,p) = exp(-x) - p # to solve p = exp(-x)
f (generic function with 1 method)

julia> fx = ZeroProblem(f, 1)
ZeroProblem{typeof(f), Int64}(f, 1)

julia> solve(fx, 1/2)  # log(2)
0.6931471805599453
```

This would be recommended, as there is no recompilation due to the function changing.

The argument `verbose=true` for `init` instructs that steps to be logged;

The iterator interface allows for the creation of hybrid solutions, for example, this is essentially how `Order0` is constructed (`Order0` follows secant steps until a bracket is identified, after which it switches to a bracketing algorithm.)

```jldoctest find_zero
julia> function order0(f, x)
           fx = ZeroProblem(f, x)
           p = init(fx, Roots.Secant())
           xᵢ,st = ϕ = iterate(p)
           while ϕ !== nothing
               xᵢ, st = ϕ
               state, ctr = st
               fᵢ₋₁, fᵢ = state.fxn0, state.fxn1
               if sign(fᵢ₋₁)*sign(fᵢ) < 0 # check for bracket
                   x0 = (state.xn0, state.xn1)
                   fx′ = ZeroProblem(f, x0)
                   p = init(fx′, Bisection())
                   xᵢ = solve!(p)
                   break
               end
               ϕ = iterate(p, st)
           end
           xᵢ
       end
order0 (generic function with 1 method)

julia> order0(sin, 3)
3.141592653589793
```

"""
function solve!(P::ZeroProblemIterator; verbose=false)
    M, F, state, options, l = P.M, P.F, P.state, P.options, P.logger

    val, stopped = :not_converged, false
    ctr = 1
    log_step(l, M, state; init=true)

    while !stopped
        val, stopped = assess_convergence(M, state, options)
        stopped && break
        ctr > options.maxevals && break

        state, stopped = update_state(M, F, state, options, l)

        log_step(l, M, state)
        ctr += 1
    end

    val, stopped = assess_convergence(M, state, options) # udpate val flag
    α = decide_convergence(M, F, state, options, val)

    log_convergence(l, val)
    log_method(l, M)
    log_last(l, α)
    verbose && display(l)

    α

end

# thread verbose through
"""
    solve(𝐙::ZeroProblem, args...; verbose=false, kwargs...)

Disptaches to `solve!(init(𝐙, args...; kwargs...))`. See [`solve!`](@ref) for details.
"""
CommonSolve.solve(F::ZeroProblem, args...; verbose=false, kwargs...) =
    solve!(init(F, args...; verbose=verbose, kwargs...); verbose=verbose)

## --------------------------------------------------
## Framework for setting up an iterative problem for finding a zero
#
# In McNamee & Pan (DOI:10.1016/j.camwa.2011.11.015 there are a number of
# results on efficiencies of a solution, (1/d) log_10(q)
# Those implemented here are:
# quadratic cut (Muller) .265 (in a42)
# Newton() newton = .1505   or 1/2log(2)
# Order1() secant method .20 (1/1 * log(1.6)
# FalsePostion(12) Anderson-Bjork [.226, .233]
# FalsePostion(3) (King?) .264
# A42() 0.191 but convergence guaranteed
# Order8() 8th order 4 steps: .225 (log10(8)/4
# Order16() 16th order 5 steps .240
# Order5(): 5th order, 4 steps. 0.1747

# A zero is found by specifying:
# the method to use <: AbstractUnivariateZeroMethod
# the function(s)
# the initial state through a value for x either x, [a,b], or (a,b) <: AbstractUnivariateZeroState
# the options (e.g., tolerances) <: UnivariateZeroOptions

# The minimal amount needed to add a method, is to define a Method and an update_state method.
