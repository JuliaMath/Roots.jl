## Framework for setting up an iterative problem for finding a zero
## TODO
## * a graphic of trace when verbose=true?


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
# the function(s) <: CallableFunction
# the initial state through a value for x either x, [a,b], or (a,b) <: AbstractUnivariateZeroState
# the options (e.g., tolerances) <: UnivariateZeroOptions

# The minimal amount needed to add a method, is to define a Method and an update_state method.

### Methods
abstract type AbstractUnivariateZeroMethod end
abstract type AbstractBracketing <: AbstractUnivariateZeroMethod end
abstract type AbstractNonBracketing <: AbstractUnivariateZeroMethod end
abstract type AbstractSecant <: AbstractNonBracketing end

# indicate if we expect f() to return one or multiple values (e.g. Newton)
fn_argout(::AbstractUnivariateZeroMethod) = 1


### States
abstract type  AbstractUnivariateZeroState{T,S} end

struct UnivariateZeroState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    fxn1::S
    fxn0::S
end

# init_state(M, F, state) -- convert
# basic idea to convert from N to M:
# Fₘ = copy(M, fₙ)
# stateₘ = init_state(M, Fₘ, stateₙ)
function init_state(M::AbstractUnivariateZeroState, F, state::AbstractUnivariateZeroState)
    init_state(M, F, state.xn0, state.xn1, state.fxn0, state.fxn1)
end


# init_state(M,F,x) --> call init_state(M,F,x₀,x₁,fx₀, fx₁)
function init_state(M::AbstractUnivariateZeroState, F, x)
    error("no default method")
end

# initialize from xs, fxs
function init_state(::AbstractUnivariateZeroState, F, x₀, x₁, fx₀, fx₁)
    error("no default method")
end


# init_state(M, F, x; kwargs...)
# init_state(M, xs..., fs...; kwargs...);
#
# A state holds at a minimum
#
# * the values xₙ₋₁, xₙ and f(xₙ₋₁), f(xₙ) along with
# * some method-specific values
#
# A state is initialized with `init_state(M, F, x)` which sets up xₙ₋₁, xₙ, f(xₙ₋₁), f(xₙ)
# which then calls `init_state(M, F, xₙ₋₁, xₙ, f(xₙ₋₁), f(xₙ))` to finish the initialization
# to change to a new state `init_state(M,F, state)`
#

# how many function evaluations in init_state
initial_fncalls(M::AbstractUnivariateZeroState) = @warn "initial_fncalls fix $M"

## -----
### Options
struct UnivariateZeroOptions{Q,R,S,T}
    xabstol::Q
    xreltol::R
    abstol::S
    reltol::T
    maxevals::Int
    maxfnevals::Int
    strict::Bool
end

"""
    default_tolerances(M::AbstractUnivariateZeroMethod, [T], [S])

The default tolerances for most methods are `xatol=eps(T)`,
`xrtol=eps(T)`, `atol=4eps(S)`, and `rtol=4eps(S)`, with the proper
units (absolute tolerances have the units of `x` and `f(x)`; relative
tolerances are unitless). For `Complex{T}` values, `T` is used.

The number of iterations is limited by `maxevals=40`, the number of
function evaluations is not capped.

"""
default_tolerances(M::AbstractUnivariateZeroMethod) = default_tolerances(M, Float64, Float64)
function default_tolerances(::AbstractUnivariateZeroMethod, ::Type{T}, ::Type{S}) where {T, S}
    xatol = eps(real(T)) * oneunit(real(T))
    xrtol = eps(real(T))  # unitless
    atol = 4 * eps(real(float(S))) * oneunit(real(S))
    rtol = 4 * eps(real(float(S))) * one(real(S))
    maxevals = 40
    maxfnevals = typemax(Int)
    strict = false
    (xatol, xrtol, atol, rtol, maxevals, maxfnevals, strict)
end

init_options(M::AbstractUnivariateZeroMethod,
             state::AbstractUnivariateZeroState{T,S};
             kwargs...
             ) where {T, S} = init_options(M, T, S; kwargs...)

function init_options(M, T=Float64, S=Float64; kwargs...)
    d = kwargs

    defs = default_tolerances(M, T, S)
    options = UnivariateZeroOptions(get(d, :xatol, get(d, :xabstol, defs[1])),
                                    get(d, :xrtol, get(d, :xreltol, defs[2])),
                                    get(d, :atol,  get(d, :abstol,  defs[3])),
                                    get(d, :rtol,  get(d, :reltol,  defs[4])),
                                    get(d, :maxevals,   get(d, :maxsteps, defs[5])),
                                    get(d, :maxfnevals, defs[6]),
                                    get(d, :strict,     defs[7]))
    options
end

# # reset options to default values
@deprecate init_options!(options, M) init_options(M)

## Tracks (for logging actual steps)
## when no logging this should get optimized out to avoid a branch
abstract type AbstractTracks end
struct NullTracks <: AbstractTracks end
# api
log_step(s::NullTracks, M, x, init=false) = nothing
log_steps(::NullTracks, n=1) = nothing
incfn(::NullTracks, i=1) = nothing
log_message(::NullTracks, msg) = nothing
log_convergence(::NullTracks, msg) = nothing
log_state(::NullTracks, state) = nothing
log_method(::NullTracks, method) = nothing
log_nmethod(::NullTracks, method) = nothing

# a tracks object to record tracks
mutable struct Tracks{T,S} <: AbstractTracks
    xs::Vector{T}
    fs::Vector{S}
    steps::Int
    fncalls::Int
    convergence_flag::String
    message::String
    state
    method
    nmethod
end
Tracks(T,S) = Tracks(T[],S[], 0, 0, "", "", nothing, nothing,nothing)
Tracks(s::AbstractUnivariateZeroState{T,S}) where {T, S} = Tracks(T,S)
Tracks(verbose, tracks, state) = (verbose && isa(tracks, NullTracks)) ? Tracks(state) : tracks

function log_step(l::Tracks, M::Any, state, init=nothing)
    if init != nothing
        x₀, fx₀ = state.xn0, state.fxn0
        push!(l.xs, x₀)
        push!(l.fs, fx₀)
    end
    x₁, fx₁ = state.xn1, state.fxn1
    push!(l.xs, x₁)
    push!(l.fs, fx₁)

    init != nothing && log_steps(l, 1)
    nothing
end

function log_steps(l::Tracks, n=1)
    l.steps += n
    nothing
end

function incfn(l::Tracks, i=1)
    l.fncalls += i
    nothing
end

function log_message(l::Tracks, msg)
    l.message *= msg
    nothing
end

function log_convergence(l::Tracks, msg)
    l.convergence_flag = string(msg)
    nothing
end

function log_state(l::Tracks, state)
    l.state = state
    nothing
end

function log_method(l::Tracks, method)
    l.method = method
    nothing
end
function log_nmethod(l::Tracks, method)
    l.nmethod = method
    nothing
end

Base.show(io::IO, l::Tracks) = show_trace(io, l.method, l.nmethod, l.state, l)

function show_tracks(io::IO, s::Tracks, M::AbstractUnivariateZeroMethod)
    for (i, (xi, fxi)) in enumerate(zip(s.xs, s.fs))
        println(io,@sprintf("%s = % 18.16f,\t %s = % 18.16f", "x_$(i-1)", float(xi), "fx_$(i-1)", float(fxi)))
    end
    println(io,"")
end

function show_trace(io::IO, method, N, state, tracks)

    converged = !isnan(state.xn1)
    println(io, "Results of univariate zero finding:\n")
    if converged
        println(io, "* Converged to: $(state.xn1)")
        if N == nothing || isa(method, AbstractBracketing)
            println(io, "* Algorithm: $(method)")
        else
            println(io, "* Algorithm: $(method), with possible bracketing with $N")
        end
        println(io, "* iterations: $(tracks.steps)")
        println("* function evaluations ≈ $(tracks.fncalls)")

        tracks.convergence_flag == :x_converged && println(io, "* stopped as x_n ≈ x_{n-1} using atol=xatol, rtol=xrtol")
        tracks.convergence_flag == :f_converged && state.message == "" && println(io, "* stopped as |f(x_n)| ≤ max(δ, max(1,|x|)⋅ϵ) using δ = atol, ϵ = rtol")
        tracks.message != "" && println(io, "* Note: $(tracks.message)")
    else
        println(io, "* Convergence failed: $(tracks.message)")
        println(io, "* Algorithm $(method)")
    end
    println(io, "")
    println(io, "Trace:")
    show_tracks(io, tracks, method)
end



### Functions
# A hacky means to call a function so that parameters can be passed as desired
# and the correct number of outputs are computed
struct Callable_Function{Single, Tup, F, P}
    f::F
    p::P
    function Callable_Function(M, f, p=nothing)
        Single = Val{fn_argout(M)}
        Tup = Val{isa(f, Tuple)}
        F = typeof(f)
        P = typeof(p)
        new{Single, Tup, F, P}(f, p)
    end
end
function CallableFunction(M,F::Callable_Function, p=F.p)
    CallableFunction(M, F.f, p)
end


# Callable_Function(output_arity, input_arity, F, p)
(F::Callable_Function{S,T,𝑭,P})(x) where {S <: Val{1}, T <: Val{false}, 𝑭, P<:Nothing} =
    first(F.f(x))
(F::Callable_Function{S,T,𝑭,P})(x) where {S <: Val{1}, T <: Val{false}, 𝑭, P} =
    first(F.f(x, F.p))

(F::Callable_Function{S,T,𝑭,P})(x) where {S <: Val{1}, T <: Val{true}, 𝑭, P<:Nothing} =
    first(F.f)(x)
(F::Callable_Function{S,T,𝑭,P})(x) where {S <: Val{1}, T <: Val{true}, 𝑭, P} =
    first(F.f)(x, F.p)

(F::Callable_Function{S,T,𝑭,P})(x) where {N, S <: Val{N}, T <: Val{false}, 𝑭, P<:Nothing} =
    F.f(x)[1:N]
(F::Callable_Function{S,T,𝑭,P})(x) where {N, S <: Val{N}, T <: Val{false}, 𝑭, P} =
    F.f(x, F.p)[1:N]

(F::Callable_Function{S,T,𝑭,P})(x) where {S <: Val{2}, T <: Val{true}, 𝑭, P<:Nothing} = begin
    f, f′ = (F.f[1])(x), (F.f[2])(x)
    (f, f/f′)
end
(F::Callable_Function{S,T,𝑭,P})(x) where {S <: Val{2}, T <: Val{true}, 𝑭, P} = begin
    f, f′ = (F.f[1])(x, F.p), (F.f[2])(x, F.p)
    (f, f/f′)
end

(F::Callable_Function{S,T,𝑭,P})(x) where {S <: Val{3}, T <: Val{true}, 𝑭, P<:Nothing} = begin
    f, f′, f′′ = (F.f[1])(x), (F.f[2])(x), (F.f[3])(x)
    (f, f/f′, f′/f′′)
end
(F::Callable_Function{S,T,𝑭,P})(x) where {S <: Val{3}, T <: Val{true}, 𝑭, P} = begin
    f, f′ = (F.f[1])(x, F.p), (F.f[2])(x, F.p), (F.f[3])(x, F.p)
    (f, f/f′, f′/f′′)
end

_apply(f,x) = f(x)
_apply(f,x, p) = f(x, p)


(F::Callable_Function{S,T,𝑭,P})(x) where {𝐍, S <: Val{𝐍}, T <: Val{true}, 𝑭, P<:Nothing} = begin
    fs = _apply.(F.f, Ref(x))
    Tuple(iszero(i) ? fs[1] : fs[i]/fs[i+1] for i ∈ 0:length(fs)-1)
end

(F::Callable_Function{S,T,𝑭,P})(x) where {𝐍, S <: Val{𝐍}, T <: Val{true}, 𝑭, P} = begin
    fs = _apply.(F.f, Ref(x), Ref(p))
    Tuple(iszero(i) ? fs[1] : fs[i]/fs[i+1] for i ∈ 0:length(fs)-1)
end



## Assess convergence
@inline function _is_f_approx_0(fa, a, atol, rtol, relaxed::Any)
    aa, afa = abs(a), abs(fa)
    tol = max(_unitless(atol), _unitless(aa) * rtol)
    tol = cbrt(abs(_unitless(tol)))  # relax test
    afa <= tol * oneunit(afa)
end
@inline function _is_f_approx_0(fa, a, atol, rtol)
    aa, afa = abs(a), abs(fa)
    tol = max(_unitless(atol), _unitless(aa) * rtol)
    afa <= tol * oneunit(afa)
end

"""
   Roots.assess_convergence(method, state, options)

Assess if algorithm has converged.

Return a convergence flag and a boolean indicating if algorithm has terminated (converged or not converged)

If alogrithm hasn't converged returns `false`.

If algorithm has stopped or converged, return `true` and sets one of `state.stopped`, `state.x_converged`,  `state.f_converged`, or `state.convergence_failed`; as well, a message may be set.

* `:x_converged` if `abs(xn1 - xn0) < max(xatol, max(abs(xn1), abs(xn0)) * xrtol)`

* `:f_converged` if  `|f(xn1)| < max(atol, |xn1|*rtol)`

* `:nan`, `:inf` if xn1 or fxn1 is `NaN` or an infinity

* `:not_converged` if algorithm should continue

Does not check number of steps taken nor number of function evaluations

In `find_zero`, stopped values (and x_converged) are checked for convergence with a relaxed tolerance.


"""
function assess_convergence(method::Any, state::AbstractUnivariateZeroState{T,S}, options) where {T,S}

    # return convergence_flag, boolean

    xn0, xn1 = state.xn0, state.xn1
    fxn1 = state.fxn1

    if isnan(xn1) || isnan(fxn1)
        return (:nan, true)
    end

    if isinf(xn1) || isinf(fxn1)
        return (:inf, true)
    end

    # f(xstar) ≈ xstar * f'(xstar)*eps(), so we pass in lambda
    if _is_f_approx_0(fxn1, xn1, options.abstol, options.reltol)
        return (:f_converged, true)
    end

    # stop when xn1 ~ xn.
    # in find_zeros there is a check that f could be a zero with a relaxed tolerance
    Δ = abs(xn1 - xn0)
    δ = max(options.xabstol, max(abs(xn1), abs(xn0)) * options.xreltol)
    if Δ ≤ δ
        return (:x_converged, true)
    end

    return (:not_converged, false)
end


# state has stopped, this identifies if it has converged
"""
    decice_convergence(M,F,state,options, convergence_flag)

When the algorithm terminates, this function decides the stopped value or returns NaN
"""
function decide_convergence(M::AbstractUnivariateZeroMethod,  F, state::AbstractUnivariateZeroState{T,S}, options, val) where {T,S}

    xn0, xn1 = state.xn0, state.xn1
    fxn1 = state.fxn1

    val ∈ (:f_converged, :exact_zero, :converged) && return xn1

    ## stopping is a heuristic, x_converged can mask issues
    ## if strict == false, this will also check f(xn) ~ - with a relaxed
    ## tolerance
    if options.strict
        val == :x_converged && return xn1
        _is_f_approx_0(fxn1, xn1, options.abstol, options.reltol) && return xn1
    else
        if val ∈ (:x_converged, :not_converged)
            # |fxn| <= relaxed
            _is_f_approx_0(fxn1, xn1, options.abstol, options.reltol, true) && return xn1
        end
    end


    # if (state.stopped || state.x_converged) && !(state.f_converged)
    #     ## stopped is a heuristic, x_converged can mask issues
    #     ## if strict == false, this will also check f(xn) ~ - with a relaxed
    #     ## tolerance

    #     ## are we at a crossing values?
    #     ## seems worth a check for 2 fn evals.
    #     if T <: Real && S <: Real
    #         for u in (prevfloat(xn1), nextfloat(xn1))
    #             fu = first(F(u))
    #             incfn(state)
    #             if iszero(fu) || _unitless(fu * fxn1) < 0
    #                 state.message *= "Change of sign at xn identified. "
    #                 state.f_converged = true
    #             end
    #         end
    #     end

    #     δ = maximum(_unitless, (options.abstol, options.reltol))
    #     if options.strict || iszero(δ)
    #         if state.x_converged
    #             state.f_converged = true
    #         else
    #             state.convergence_failed = true
    #         end

    #     else
    #         xstar, fxstar = state.xn1, state.fxn1
    #         if _is_f_approx_0(fxstar, xstar, options.abstol, options.reltol, :relaxed)
    #             state.xstar, state.fxstar = xstar, fxstar
    #             msg = "Algorithm stopped early, but |f(xn)| < ϵ^(1/3), where ϵ depends on xn, rtol, and atol. "
    #             state.message = state.message == "" ? msg : state.message * "\n\t" * msg
    #             state.f_converged = true
    #         else
    #             state.convergence_failed = true
    #         end
    #     end
    # end

    # if !state.f_converged
    #     state.xstar, state.fxstar = NaN*xn1, NaN*fxn1
    # end

    NaN * xn1

end

## ----

"""

    find_zero(fs, x0, M, [N::AbstractBracketing]; kwargs...)

Interface to one of several methods for find zeros of a univariate function.

# Initial starting value

For most methods, `x0` is a scalar value indicating the initial value
in the iterative procedure. (Secant methods can have a tuple specify
their initial values.) Values must be a subtype of `Number` and have
methods for `float`, `real`, and `oneunit` defined.

For bracketing intervals, `x0` is specified as a tuple or a vector. A bracketing interval, (a,b), is one where f(a) and f(b) have different signs.

# Specifying a method

A method is specified to indicate which algorithm to employ:

* There are methods for bisection where a bracket is specified: [`Bisection`](@ref), [`Roots.A42`](@ref), [`Roots.AlefeldPotraShi`](@ref), [`FalsePosition`](@ref)

* There are several derivative-free methods: cf. [`Order0`](@ref), [`Order1`](@ref) (secant method), [`Order2`](@ref) ([`Roots.Steffensen`](@ref)), [`Order5`](@ref), [`Order8`](@ref), and [`Order16`](@ref), where the number indicates the order of the convergence. Methods [`Roots.Order1B`](@ref) and [`Roots.Order2B`](@ref) implement methods useful when the desired zero has a multiplicity.

* There are some classical methods where derivatives are required: [`Roots.Newton`](@ref), [`Roots.Halley`](@ref), [`Roots.Schroder`](@ref).

* The family [`Roots.LithBoonkkampIJzerman{S,D}`](@ref) for different `S` and `D` uses a linear multistep method root finder. The (2,0) method is the secant method, (1,1) is Newton's methods.

For more detail, see the help page for each method (e.g., `?Order1`). Many methods are not exported, so much be qualified with module name, as in `?Roots.Schroder`.

If no method is specified, the default method depends on `x0`:

* If `x0` is a scalar, the default is the slower, but more robust `Order0` method.

* If `x0` is a tuple or vector indicating a *bracketing* interval, the `Bisection` method is used. (The exact algorithm depends on the number type, the tolerances, and `verbose`.)

# Specifying the function

The function(s) are passed as the first argument.

For the few methods that use one or more derivatives (`Newton`, `Halley`,
`Schroder`, `LithBoonkkampIJzerman(S,D)`, and  `Order5Derivative`) a
tuple of functions is used.

# Optional arguments (tolerances, limit evaluations, tracing)

* `xatol` - absolute tolerance for `x` values. Passed to `isapprox(x_n, x_{n-1})`
* `xrtol` - relative tolerance for `x` values. Passed to `isapprox(x_n, x_{n-1})`
* `atol`  - absolute tolerance for `f(x)` values.
* `rtol`  - relative tolerance for `f(x)` values.
* `maxevals`   - limit on maximum number of iterations
* `strict` - if `false` (the default), when the algorithm stops, possible zeros are checked with a relaxed tolerance
* `verbose` - if `true` a trace of the algorithm will be shown on successful completion. See the internal `Tracks` object to save this trace.

See the help string for `Roots.assess_convergence` for details on
convergence. See the help page for `Roots.default_tolerances(method)`
for details on the default tolerances.

In general, with floating point numbers, convergence must be
understood as not an absolute statement. Even if mathematically α is
an answer and xstar the floating point realization, it may be that
`f(xstar) - f(α)  ≈ xstar ⋅  f'(α) ⋅ eps(α)`, so tolerances must be
appreciated, and at times specified.

For the `Bisection` methods, convergence is guaranteed, so the tolerances are set to be 0 by default.

If a bracketing method is passed in after the method specification, if
a bracket is identified during the algorithm, the bracketing method
will then be used to identify the zero. This is what `Order0` does by
default, with initials secant steps and then `AlefeldPotraShi` steps
should a bracket be encountered.

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
3.1415926535897936
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
ERROR: Roots.ConvergenceFailed("Stopped at: xn = 2.991488255523429. Increment `Δx` has issues. Too many steps taken. ")
[...]

julia> fn = x -> (sin(x)*cos(x) - x^3 + 1)^9;

julia> x0, xstar = 1.0,  1.112243913023029;

julia> find_zero(fn, x0, Order2()) ≈ xstar
true

julia> find_zero(fn, x0, Order2(), maxevals=3)    # Roots.ConvergenceFailed: 26 iterations needed, not 3
ERROR: Roots.ConvergenceFailed("Stopped at: xn = 1.0482748172022405. Too many steps taken. ")
```

Passing `verbose=true` will show details on the steps of the algorithm:

```jldoctest find_zero
julia> find_zero(x->sin(x), 3.0, Order2(), verbose=true)   # 2 iterations
Results of univariate zero finding:

* Converged to: 3.1415926535897936
* Algorithm: Order2()
* iterations: 2
* function evaluations: 5
* stopped as |f(x_n)| ≤ max(δ, max(1,|x|)⋅ϵ) using δ = atol, ϵ = rtol

Trace:
x_0 =  3.0000000000000000,	 fx_0 =  0.1411200080598672
x_1 =  3.1425464815525403,	 fx_1 = -0.0009538278181169
x_2 =  3.1415926535897936,	 fx_2 = -0.0000000000000003

3.1415926535897936
```

!!! note
    See [`solve`](@ref) and [`ZeroProblem`](@ref) for an alternate interface.
"""
function find_zero(fs, x0, M::AbstractUnivariateZeroMethod;
                   p = nothing,
                   verbose=false,
                   tracks::AbstractTracks=NullTracks(),
                   kwargs...)

    Z = ZeroProblem(fs, x0)
    ZPI = init(Z, M, p;
               verbose=verbose, tracks=tracks,
               kwargs...)

    xstar = solve!(ZPI; verbose=verbose)

    isnan(xstar) && throw(ConvergenceFailed("Stopped"))

    xstar

end

# defaults when method is not specified
# if a number, use Order0
# O/w use a bracketing method of an assumed iterable
find_zero(f, x0::T; kwargs...)  where {T <: Number} = find_zero(f, x0, Order0(); kwargs...)
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

# Examples

```
f(x) = sin(x)
xs = (3.0, 4.0)
fxs = f.(xs)
if prod((sign∘f).(xs)) < 0 # check bracket
    M = Bisection()
    state = Roots._init_state(M, f, xs, fxs) # reuse fxs from test
    find_zero(M, f, state)
end
```

!!! note
    To be deprecated
"""
function find_zero(M::AbstractUnivariateZeroMethod,
                   F,
                   state::AbstractUnivariateZeroState,
                   options::UnivariateZeroOptions=init_options(M, state),
                   l::AbstractTracks=NullTracks()
                   )
    solve!(init(M, F, state, options, l))
end



## ---------------

## Create an Iterator interface
# returns NaN, not an error, if there are issues

"""
    ZeroProblem{F,X}

A container for a function and initial guess passed to an iterator to be solved by `find_zero!` or `solve!`.
"""
struct ZeroProblem{F,X}
    F::F
    x₀::X
end



## The actual iterating object
struct ZeroProblemIterator{M,F,S,O,L}
    M::M
    F::F
    state::S
    options::O
    logger::L
end

## Initialize a Zero Problem Iterator
## init(Z,p)
## init(Z,M,p)
## init(M,F,state, [options], [logger])
## want p to be positional, not named⁺
function init(𝑭𝑿::ZeroProblem, M::AbstractUnivariateZeroMethod, p′ = nothing;
              p = nothing,
              verbose::Bool=false,
              tracks = NullTracks(),
              kwargs...)

    F = Callable_Function(M, 𝑭𝑿.F, p === nothing ? p′ : p)  #⁺
    state = init_state(M, F, 𝑭𝑿.x₀)
    options = init_options(M, state; kwargs...)
    l = Tracks(verbose, tracks, state)
    incfn(l, initial_fncalls(M))
    ZeroProblemIterator(M,F,state,options,l)

end

function init(𝑭𝑿::ZeroProblem, p′=nothing; kwargs...)
    M = length(𝑭𝑿.x₀) == 1 ? Secant() : Bisection()
    init(𝑭𝑿, M, p′; kwargs...)
end

function init(M::AbstractUnivariateZeroMethod, F,
              state::AbstractUnivariateZeroState,
              options::UnivariateZeroOptions=init_options(M, state),
              l::AbstractTracks=NullTracks())
    ZeroProblemIterator(M, Callable_Function(M,F), state, options, l)
end

# Optional iteration interface to handle looping
# errors on non-convergence, unlike solve!(init(...)) which returns NaN
# allocates...
# This should be deprecated
function Base.iterate(P::ZeroProblemIterator, st=nothing)
    ## st = (val, (state, ctr, val))
    if st == nothing
        state = P.state
        ctr = 1
    else
        (state, ctr, val) = st
        ctr += 1
    end
    M, options, l = P.M, P.options, P.logger
    val,stopped = :no_convergence, false
    while true
        val, stopped = assess_convergence(M, state, options)
        stopped && break
        if ctr > options.maxevals
            stopped = true
            break
        end
        state, stopped = update_state(M, P.F, state, options, l)
        log_step(l, M, state)
        break
    end

    if stopped
        xstar = decide_convergence(P.M, P.F, state, P.options, val)
        isnan(xstar) && throw(ConvergenceFailed("Stopped at $(state.xn1)"))
        nothing
    else
        return (state.xn1, (state, ctr, val))
    end
end

#log_step(P::ZeroProblemIterator) = log_step(P.logger, P.M, P.state)
#decide_convergence(P::ZeroProblemIterator) =  decide_convergence(P.M, P.F, P.state, P.options)

#"Get last element in the iteration, which is xstar, or throw a warning"
#function Base.last(P::ZeroProblemIterator)
#    state.convergence_failed && @warn "The problem failed to converge"
    #    !(state.x_converged || state.f_converged) && @warn "The problem has not converged. Try calling `solve! first"
#    P.state.xstar
#end

# """
#     tracks(P::ZeroProblemIterator)

# Show trace of output when `verbose=true` is specified to the problem
# """
# function tracks(P::ZeroProblemIterator{M,F,S,O,L}) where {M,F,S,O,L<:Tracks}
#      show_trace(P.M, nothing, P.state, P.logger)
# end
# tracks(P::ZeroProblemIterator) = error("Set verbose=true when specifying the problem to see the tracks")
# function show_trace(P::ZeroProblemIterator{M,F,S,O,L}) where {M,F,S,O,L<:Tracks}
#      show_trace(P.M, nothing, P.state, P.logger)
# end

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
* `solve!` to iterate to convergence. (Also [`find_zero!`](@ref).)

Returns `NaN`, not an error, when the problem can not be solved.

## Examples:

```
fx = ZeroProblem(sin, 3)
solve(fx)
```

Or, if the iterator is required

```
fx = ZeroProblem(sin, 3)
problem = init(fx)
solve!(fx)
```

The default method is `Order1()`, when  `x0` is a number, or `Bisection()` when `x0` is an iterable with 2 or more values.


A second position argument for `solve` or `init` is used to specify a different method; keyword arguments can be used to adjust the default tolerances.


```
fx = ZeroProblem(sin,3)
solve(fx, Order5(), atol=1/100)
```

The above is equivalent to:

```
fx = ZeroProblem(sin, 3)
problem = init(fx, Order5(), atol=1/100)
solve!(problem)
```

The  argument `p` may be used if the function(s) to be solved depend on a parameter in their second positional argument (e.g., `f(x,p)`). For example

```
f(x,p) = exp(-x) - p # to solve p = exp(-x)
fx = ZeroProblem(f, 1)
solve(fx, 1/2)  # log(2)
```

This would be recommended, as there is no recompilation due to the function changing.

The argument `verbose=true` for `init` instructs that steps to be logged;

The iterator interface allows for the creation of hybrid solutions, for example, this is essentially how `Order0` is constructed (`Order0` follows secant steps until a bracket is identified, after which is switches to a bracketing algorithm.)

```
function order0(f, x)
    fx = ZeroProblem(f, x)
    p = init(fx, Roots.Secant())
    xᵢ,st = ϕ = iterate(p)
    while ϕ != nothing
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
```

"""
function solve!(P::ZeroProblemIterator; verbose=false)
    M, F, state, options, l = P.M, P.F, P.state, P.options, P.logger
    val, stopped = :not_converged, false
    ctr = 1
    log_step(l, M, state, :init)
    while !stopped
        val, stopped = assess_convergence(M, state, options)
        stopped && break
        ctr > options.maxevals && break
        #ctr * fn_argout(M) > options.maxfnevals && break
        state, stopped = update_state(M, F, state, options, l)
        log_step(l, M, state)
        ctr += 1
    end

    if !isa(l, NullTracks)
        log_convergence(l, val); log_state(l, state); log_method(l, M)
        verbose && display(l)
    end

    val, stopped = assess_convergence(M, state, options) # udpate val flag
    decide_convergence(M, F, state, options, val)

end


## -----
## deprecate this interface at some time.
@deprecate find_zero!(P::ZeroProblemIterator) solve!(P)
# """
#     find_zero!(P::ZeroProblemIterator)

# An alternate interface to `find_zero` whereby a problem is created with `ZeroProblemIterator` and solved
# with `find_zero!`. The generic [`solve!`](@ref) method is recommened for familiarity.


# ```jldoctest find_zero
# julia> using Roots

# julia> P = ZeroProblem(Order1(), sin, 3, verbose=true);

# julia> find_zero!(P)
# 3.141592653589793

# julia> last(P)
# 3.141592653589793

# julia> Roots.tracks(P) # possible when `verbose=true` is specified
# Results of univariate zero finding:

# * Converged to: 3.141592653589793
# * Algorithm: Roots.Secant()
# * iterations: 4
# * function evaluations: 6
# * stopped as |f(x_n)| ≤ max(δ, max(1,|x|)⋅ϵ) using δ = atol, ϵ = rtol

# Trace:
# x_0 =  3.0000000000000000,	 fx_0 =  0.1411200080598672
# x_1 =  3.1425464815525403,	 fx_1 = -0.0009538278181169
# x_2 =  3.1415894805773834,	 fx_2 =  0.0000031730124098
# x_3 =  3.1415926535902727,	 fx_3 = -0.0000000000004795
# x_4 =  3.1415926535897931,	 fx_4 =  0.0000000000000001
# ```
# """
#function find_zero!(P::ZeroProblemIterator)
#    solve!(P)
#end


@deprecate ZeroProblem(M::AbstractUnivariateZeroMethod,fs,x0;kwargs...) init(ZeroProblem(fs, x0), M;kwargs...)

# """
#     ZeroProblem(M, fs, x0; verbose=false, kwargs...)

# Setup an interator interface for the zero problem. Call `find_zero!` or solve! to solve.

# * `M`: Method. A non-hybrid method (not `Order0`).
# * `fs`: a function or for some methods a tuple of functions
# * `x0`: the initial guess
# * `verbose`: if `true`, then calling `Roots.tracks` on the output will show the steps on the algorithm
# * `kwargs`: passed to `Roots.init_options` to adjust tolerances

# """
# function ZeroProblem(M::AbstractUnivariateZeroMethod,
#                      fs,
#                      x0;
#                      verbose=false,
#                      kwargs...)
#     fx = ZeroProblem(fs, x0)
#     problem = init(fx, M; verbose=verbose, kwargs...)

# end
