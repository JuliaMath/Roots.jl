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
abstract type  AbstractUnivariateZeroState end

mutable struct UnivariateZeroState{T,S} <: AbstractUnivariateZeroState where {T,S}
    xn1::T
    xn0::T
    xstar::T
    m::Vector{T}
    fxn1::S
    fxn0::S
    fxstar::S
    fm::Vector{S}
    steps::Int
    fnevals::Int
    stopped::Bool             # stopped, butmay not have converged
    x_converged::Bool         # converged via |x_n - x_{n-1}| < ϵ
    f_converged::Bool         # converged via |f(x_n)| < ϵ
    convergence_failed::Bool
    message::String
end

incfn(o::AbstractUnivariateZeroState, k=1)    = o.fnevals += k
incsteps(o::AbstractUnivariateZeroState, k=1) = o.steps += k

# initialize state for most methods
function init_state(method::Any, fs, x)

    x1 = float(x)
    fx1 = fs(x1); fnevals = 1
    T, S = eltype(x1), eltype(fx1)
    zT, zS =  oneunit(x1) * (0*x1)/(0*x1), oneunit(fx1) * (0*fx1)/(0*fx1)
    state = UnivariateZeroState(x1, zT, zT/Zt*oneunit(x1), T[],
                                fx1, zS, zS, S[],
                                0, fnevals,
                                false, false, false, false,
                                "")
    state
end

## This function is used to reset the state to an initial value
## As initializing a state is somewhat costly, this can be useful when many
## function calls would be used.
## An example usage, might be:
# M = Order1()
# state = Roots.init_state(M, f1, 0.9)
# options = Roots.init_options(M, state)
# out = zeros(Float64, n)
# for (i,x0) in enumerate(linspace(0.9, 1.0, n))
#    Roots.init_state!(state, M, f1, x0)
#    out[i] = find_zero(M, f1, state, options)
# end
# init_state has a method call variant
function init_state!(state::UnivariateZeroState{T,S}, x1::T, x0::T, m::Vector{T}, y1::S, y0::S, fm::Vector{S}) where {T,S}
    state.xn1 = x1
    state.xn0 = x0
    empty!(state.m); append!(state.m,  m)
    state.m = m
    state.fxn1 = y1
    state.fxn0 = y0
    empty!(state.fm); append!(state.fm,fm)
    state.steps = 0
    state.fnevals = 2
    state.stopped = false
    state.x_converged = false
    state.f_converged = false
    state.convergence_failed = false
    state.message = ""

    nothing
end

function init_state!(state::UnivariateZeroState{T,S}, ::AbstractUnivariateZeroMethod, fs, x) where {T, S}
    x1 = float(x)
    fx1 = fs(x1)

    init_state!(state, x1, oneunit(x1) * (0*x1)/(0*x1), T[],
                fx1, oneunit(fx1) * (0*fx1)/(0*fx1), S[])
end


function Base.copy(state::UnivariateZeroState{T,S}) where {T, S}
    UnivariateZeroState(state.xn1, state.xn0, state.xstar, copy(state.m),
                        state.fxn1, state.fxn0, state.fxstar, copy(state.fm),
                        state.steps, state.fnevals,
                        state.stopped, state.x_converged,
                        state.f_converged, state.convergence_failed,
                        state.message)
end

function Base.copy!(dst::UnivariateZeroState{T,S}, src::UnivariateZeroState{T,S}) where {T,S}
    dst.xn1 = src.xn1
    dst.xn0 = src.xn0
    empty!(dst.m); append!(dst.m, src.m)
    dst.fxn1 = src.fxn1
    dst.fxn0 = src.fxn0
    empty!(dst.fm); append!(dst.fm, src.fm)
    dst.steps = src.steps
    dst.fnevals = src.fnevals
    dst.stopped = src.stopped
    dst.x_converged = src.x_converged
    dst.f_converged = src.f_converged
    dst.convergence_failed = src.convergence_failed
    dst.message = src.message
    nothing
end

### Options
mutable struct UnivariateZeroOptions{Q,R,S,T}
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
    atol = 4.0 * eps(real(float(S))) * oneunit(real(S))
    rtol = 4.0 * eps(real(float(S))) * one(real(S))
    maxevals = 40
    maxfnevals = typemax(Int)
    strict = false
    (xatol, xrtol, atol, rtol, maxevals, maxfnevals, strict)
end

function init_options(M::AbstractUnivariateZeroMethod,
                      state::UnivariateZeroState{T,S};
                      kwargs...
                      ) where {T, S}

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

# reset options to default values
function init_options!(options::UnivariateZeroOptions{Q,R,S,T}, M::AbstractUnivariateZeroMethod) where {Q, R, S, T}

    defs = default_tolerances(M, Q, S)
    options.xabstol = defs[1]
    options.xreltol = defs[2]
    options.abstol = defs[3]
    options.reltol = defs[4]
    options.maxevals = defs[5]
    options.maxfnevals = defs[6]
    options.strict = defs[7]

    nothing
end


## Tracks (for logging actual steps)
## when no logging this should get optimized out to avoid a branch
abstract type AbstractTracks end
struct NullTracks <: AbstractTracks end
# api
log_step(s::NullTracks, M, x, init=false) = nothing
log_step(::Nothing, M, x, init=false) = nothing

mutable struct Tracks{T,S} <: AbstractTracks
xs::Vector{T}
fs::Vector{S}
end
Tracks(s::UnivariateZeroState{T,S}) where {T, S} = Tracks(T[],S[])
Tracks(verbose, tracks, state) = (verbose && isa(tracks, NullTracks)) ? Tracks(state) : tracks


log_step(s::Tracks, M::Any, o, ::Any) = log_step(s, M, o)

function log_step(s::Tracks, M::Any, o)
    push!(s.xs, o.xn1)
    push!(s.fs, o.fxn1)
    nothing
end
function show_tracks(s::Tracks, M::AbstractUnivariateZeroMethod)
    for (i, (xi, fxi)) in enumerate(zip(s.xs, s.fs))
        println(@sprintf("%s = % 18.16f,\t %s = % 18.16f", "x_$(i-1)", float(xi), "fx_$(i-1)", float(fxi)))
    end
    println("")
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

If alogrithm hasn't converged returns `false`.

If algorithm has stopped or converged, return `true` and sets one of `state.stopped`, `state.x_converged`,  `state.f_converged`, or `state.convergence_failed`; as well, a message may be set.

* `state.x_converged = true` if `abs(xn1 - xn0) < max(xatol, max(abs(xn1), abs(xn0)) * xrtol)`

* `state.f_converged = true` if  `|f(xn1)| < max(atol, |xn1|*rtol)`

* `state.convergence_failed = true` if xn1 or fxn1 is `NaN` or an infinity

* `state.stopped = true` if the number of steps exceed `maxevals` or the number of function calls exceeds `maxfnevals`.

In `find_zero`, stopped values (and x_converged) are checked for convergence with a relaxed tolerance.


"""
function assess_convergence(method::Any, state::UnivariateZeroState{T,S}, options) where {T,S}

    xn0, xn1 = state.xn0, state.xn1
    fxn1 = state.fxn1

    if (state.x_converged || state.f_converged) #|| state.stopped)
        if isnan(state.xstar)
            state.xstar, state.fxstar =  xn1, fxn1
        end
        return true
    end

    if isnan(xn1) || isnan(fxn1)
        state.convergence_failed = true
        state.message *= "NaN produced by algorithm. "
        return true
    end

    if isinf(xn1) || isinf(fxn1)
        state.convergence_failed = true
        state.message *= "Inf produced by algorithm. "
        return true
    end

    # f(xstar) ≈ xstar * f'(xstar)*eps(), so we pass in lambda
    if _is_f_approx_0(fxn1, xn1, options.abstol, options.reltol)
        state.xstar, state.fxstar = xn1, fxn1
        state.f_converged = true
        return true
    end

    # stop when xn1 ~ xn.
    # in find_zeros there is a check that f could be a zero with a relaxed tolerance
    if abs(xn1 - xn0) < max(options.xabstol, max(abs(xn1), abs(xn0)) * options.xreltol)
        state.xstar, state.fxstar = xn1, fxn1
        state.message *= "x_n ≈ x_{n-1}. "
        state.x_converged = true
        return true
    end

    check_steps_fnevals(state, options) && return true

    return false
end

## too many steps or function evaluations
function check_steps_fnevals(state, options)

    if state.steps > options.maxevals
        state.stopped = true
        state.message *= "Too many steps taken. "
        return true
    end

    if state.fnevals > options.maxfnevals
        state.stopped = true
        state.message *= "Too many function evaluations taken. "
        return true
    end

    return false
end

# state has stopped, this identifies if it has converged
function decide_convergence(M::AbstractUnivariateZeroMethod,  F, state::UnivariateZeroState{T,S}, options) where {T,S}
    xn1 = state.xstar
    fxn1 = state.fxstar

    if (state.stopped || state.x_converged) && !(state.f_converged)
        ## stopped is a heuristic, x_converged can mask issues
        ## if strict == false, this will also check f(xn) ~ - with a relaxed
        ## tolerance

        ## are we at a crossing values?
        ## seems worth a check for 2 fn evals.
        if T <: Real && S <: Real
            for u in (prevfloat(xn1), nextfloat(xn1))
                fu = first(F(u))
                if iszero(fu) || _unitless(fu * fxn1) < 0
                    state.message *= "Change of sign at xn identified. "
                    state.f_converged = true
                end
            end
        end

        δ = maximum(_unitless, (options.abstol, options.reltol))
        if options.strict || iszero(δ)
            if state.x_converged
                state.f_converged = true
            else
                state.convergence_failed = true
            end

        else
            xstar, fxstar = state.xn1, state.fxn1
            if _is_f_approx_0(fxstar, xstar, options.abstol, options.reltol, :relaxed)
                state.xstar, state.fxstar = xstar, fxstar
                msg = "Algorithm stopped early, but |f(xn)| < ϵ^(1/3), where ϵ depends on xn, rtol, and atol. "
                state.message = state.message == "" ? msg : state.message * "\n\t" * msg
                state.f_converged = true
            else
                state.convergence_failed = true
            end
        end
    end

    if !state.f_converged
        state.xstar, state.fxstar = NaN*xn1, NaN*fxn1
    end

    state.xstar

end

## ----

function show_trace(method, N, state, tracks)
    converged = state.x_converged || state.f_converged
    println("Results of univariate zero finding:\n")
    if converged
        println("* Converged to: $(state.xn1)")
        if N == nothing || isa(method, AbstractBracketing)
            println("* Algorithm: $(method)")
        else
            println("* Algorithm: $(method), with possible bracketing with $N")
        end
        println("* iterations: $(state.steps)")
        println("* function evaluations: $(state.fnevals)")
        state.x_converged && println("* stopped as x_n ≈ x_{n-1} using atol=xatol, rtol=xrtol")
        state.f_converged && state.message == "" && println("* stopped as |f(x_n)| ≤ max(δ, max(1,|x|)⋅ϵ) using δ = atol, ϵ = rtol")
        state.message != "" && println("* Note: $(state.message)")
    else
        println("* Convergence failed: $(state.message)")
        println("* Algorithm $(method)")
    end
    println("")
    println("Trace:")
    show_tracks(tracks, method)
end

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
`Schroder`, `LithBoonkkampIJzerman(S,D)`, and optionally `Order5`) a
tuple of functions is used.

# Optional arguments (tolerances, limit evaluations, tracing)

* `xatol` - absolute tolerance for `x` values. Passed to `isapprox(x_n, x_{n-1})`
* `xrtol` - relative tolerance for `x` values. Passed to `isapprox(x_n, x_{n-1})`
* `atol`  - absolute tolerance for `f(x)` values.
* `rtol`  - relative tolerance for `f(x)` values.
* `maxevals`   - limit on maximum number of iterations
* `maxfnevals` - limit on maximum number of function evaluations
* `strict` - if `false` (the default), when the algorithm stops, possible zeros are checked with a relaxed tolerance
* `verbose` - if `true` a trace of the algorithm will be shown on successful completion.

See the help string for `Roots.assess_convergence` for details on
convergence. See the help page for `Roots.default_tolerances(method)`
for details on the default tolerances.

In general, with floating point numbers, convergence must be
understood as not an absolute statement. Even if mathematically α is
an answer and xstar the floating point realization, it may be that
`f(xstar) - f(α)  ≈ xstar ⋅  f'(α) ⋅ eps(α)`, so tolerances must be
appreciated, and at times specified.

For the `Bisection` methods, convergence is guaranteed, so the tolerances are set to be 0 by default.

If a bracketing method is passed in after the method specification, when a bracket is found, the bracketing method will used to identify the zero. This is what `Order0` does by default, with a secant step initially and the `A42` method when a bracket is  encountered.

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
ERROR: Roots.ConvergenceFailed("Stopped at: xn = 2.9921081308371615. Too many steps taken. ")

julia> fn = x -> (sin(x)*cos(x) - x^3 + 1)^9;

julia> x0, xstar = 1.0,  1.112243913023029;

julia> find_zero(fn, x0, Order2()) ≈ xstar
true

julia> find_zero(fn, x0, Order2(), maxevals=3)    # Roots.ConvergenceFailed: 26 iterations needed, not 3
ERROR: Roots.ConvergenceFailed("Stopped at: xn = 1.0482748172022405. Too many steps taken. ")
```

Tracing output.

```jldoctest find_zero
julia> find_zero(x->sin(x), 3.0, Order2(), verbose=true)   # 3 iterations
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

```jldoctest find_zero
julia> find_zero(x->sin(x)^5, 3.0, Order2(), verbose=true) # 22 iterations
Results of univariate zero finding:

* Converged to: 3.1405349388309682
* Algorithm: Order2()
* iterations: 22
* function evaluations: 46
* stopped as |f(x_n)| ≤ max(δ, max(1,|x|)⋅ϵ) using δ = atol, ϵ = rtol

Trace:
x_0 =  3.0000000000000000,	 fx_0 =  0.0000559684091879
x_1 =  3.0285315910604353,	 fx_1 =  0.0000182783542076
x_2 =  3.0512479368872216,	 fx_2 =  0.0000059780426727
[...]
x_21 =  3.1402867582538070,	 fx_21 =  0.0000000000000038
x_22 =  3.1405349388309682,	 fx_22 =  0.0000000000000013

3.1405349388309682

julia> find_zero(x->sin(x)^5, 3.0, Roots.Order2B(), verbose=true) # 2 iterations
Results of univariate zero finding:

* Converged to: 3.1397074174874358
* Algorithm: Roots.Order2B()
* iterations: 41
* function evaluations: 85
* Note: Estimate for multiplicity had issues. Too many steps taken.
	  Algorithm stopped early, but |f(xn)| < ϵ^(1/3), where ϵ depends on xn, rtol, and atol.

Trace:
x_0 =  3.0000000000000000,	 fx_0 =  0.0000559684091879
x_1 =  3.1397074174874358,	 fx_1 =  0.0000000000000238
[...]
x_40 =  3.1397074174874358,	 fx_40 =  0.0000000000000238
x_41 =  3.1397074174874358,	 fx_41 =  0.0000000000000238

3.1397074174874358
```

!!! note
    See [`ZeroProblem`](@ref) for an interator interface to the underlying algorithms.
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

    xstar = solve!(ZPI)
    verbose && show_trace(M, nothing, ZPI.state, ZPI.logger)
    isnan(xstar) && throw(ConvergenceFailed("Stopped at: xn = $(ZPI.state.xn1). $(ZPI.state.message)"))

    return xstar

end

# defaults when method is not specified
# if a number, use Order0
# O/w use a bracketing method of an assumed iterable
find_zero(f, x0::T; kwargs...)  where {T <: Number} = find_zero(f, x0, Order0(); kwargs...)
find_zero(f, x0; kwargs...) = find_zero(f, x0, Bisection(); kwargs...)

# todo: consolidate this with find_zero(M,N,f, x0)...
function find_zero(fs, x0,
                   M::AbstractUnivariateZeroMethod,
                   N::AbstractBracketing;
                   p = nothing,
                   tracks::AbstractTracks=NullTracks(),
                   verbose=false,
                   kwargs...)

    F = Callable_Function(M, fs, p)
    state = init_state(M, F, x0)
    options = init_options(M, state; kwargs...)
    l = Tracks(verbose, tracks, state)

    xstar = find_zero(M, N, F, state, options, l)
    verbose &&  show_trace(M, N, state, l)

    isnan(xstar) && throw(ConvergenceFailed("Stopped at: xn = $(state.xn1). $(state.message)"))

    return xstar

end




## find_zero implementation
## Updates state, could be `find_zero!(state, M, F, options, l)`...
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

# Robust version using some tricks: idea from algorithm described in
# [The SOLVE button from the
# HP-34]C(http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf).
# * use bracketing method if one identifed
# * limit steps so as not too far or too near the previous one
# * if not decreasing, use a quad step upto 4 times to bounce out of trap, if possible
# First uses M, then N if bracket is identified
function find_zero(M::AbstractUnivariateZeroMethod,
                   N::AbstractBracketing,
                   F,
                   state::AbstractUnivariateZeroState,
                   options::UnivariateZeroOptions,
                   l::AbstractTracks=NullTracks()
                   )


    log_step(l, M, state, :init)
    state0 = copy(state)
    quad_ctr = 0

    while true

        if assess_convergence(M, state, options)
            break
        end

        copy!(state0, state)
        update_state(M, F, state0, options) # state0 is proposed step

        adj = false
        ## did we find a zero or a bracketing interval?
        if iszero(state0.fxn1)
            copy!(state, state0)
            state.xstar, state.fxstar = state.xn1, state.fxn1
            state.f_converged = true
            break
        elseif sign(state0.fxn0) * sign(state0.fxn1) < 0
            copy!(state, state0)
            a, b = state0.xn0, state0.xn1 # could save some fn calls here
            run_bisection(N, F, (a, b), state, options)
            break
        end


        ## did we move too far?
        r, a, b = state0.xn1, state.xn0, state.xn1
        Δr = abs(r - b)
        Δx = abs(b - a)
        ts,TB = 1e-3, 1e2 # too small, too big
        if  Δr >= TB * Δx
            adj = true
            r = b + sign(r-b) * TB * Δx  ## too big
            state0.xn1 = r
            state0.fxn1 = F(r)
            incfn(state)
        elseif Δr  <= ts *  Δx
            adj = true
            r = b + sign(r - b) * ts * Δx
            fr = F(r)
            incfn(state)
            state0.xn1 = r
            state0.fxn1 = fr
        end

        # a sign change after shortening?
        if sign(state.fxn1) * sign(state0.fxn1) < 0
            state.xn0, state.fxn0 = state.xn1, state.fxn1
            state.xn1, state.fxn1 = state0.xn1, state0.fxn1
            a, b = state.xn0, state.xn1
            run_bisection(N, F, (a, b), state, options)
            break
        end


        ## did we improve?
        if adj || abs(state0.fxn1) < abs(state.fxn1)
            if isnan(state0.xn1) || isnan(state0.fxn1) || isinf(state0.xn1) || isinf(state0.fxn1)
                break
            end
            copy!(state, state0)
            log_step(l, M, state)
            incsteps(state)
            quad_ctr = 0
            continue
        end

        ## try quad_vertex, unless that has gotten old
        if quad_ctr > 4
            copy!(state, state0)
            state.stopped = true
            break
        else
            quad_ctr += 1
            r = quad_vertex(state0.xn1, state0.fxn1, state.xn1, state.fxn1, state.xn0, state.fxn0)

            if isnan(r) || isinf(r)
                copy!(state, state0)
            else

                fr = F(r)
                incfn(state)

                state0.xn1 = r
                state0.fxn1 = fr
                incfn(state)
                copy!(state, state0)
            end
        end

        log_step(l, M, state)
        incsteps(state)
    end

    decide_convergence(M, F, state, options)
end


## --


## ---- Hybrid Method ----
# Switch to bracketing method
#M = Bisection()  # exact for floating point
#M = AlefeldPotraShi() # *usually* exact
#M = Brent()          # a bit faster, but not always convergent, as implemented (cf. RootTesting)
run_bisection(f, xs, state, options) = run_bisection(AlefeldPotraShi(), f, xs, state, options)
function run_bisection(N::AbstractBracketing, f, xs, state, options)
    steps, fnevals = state.steps, state.fnevals
    f = Callable_Function(N, f)
    init_state!(state, N, f, xs)
    state.steps += steps
    state.fnevals += fnevals # could avoid 2 fn calls, with fxs
    init_options!(options, N)
    find_zero(N, f, state, options)
    a, b = xs
    u,v = a > b ? (b, a) : (a, b)
    state.message *= "Bracketing used over ($u, $v), those steps not shown. "
    return nothing
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
              verbose=false,
              tracks = NullTracks(),
              kwargs...)

    F = Callable_Function(M, 𝑭𝑿.F, p === nothing ? p′ : p)  #⁺
    state = init_state(M, F, 𝑭𝑿.x₀)
    options = init_options(M, state; kwargs...)
    l = Tracks(verbose, tracks, state)
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

# Iteration interface to handle looping
function Base.iterate(P::ZeroProblemIterator, st=nothing)

    M, F, state, options, l = P.M, P.F, P.state, P.options, P.logger
    st == nothing && log_step(l, M, state, :init)

    if assess_convergence(M, state, options)
        decide_convergence(M,F,state,options)
        return  nothing
    end

    update_state(M, F, state, options)
    log_step(l, M, state)
    incsteps(state)
    (state.xn1, false)

end

log_step(P::ZeroProblemIterator) = log_step(P.logger, P.M, P.state)
decide_convergence(P::ZeroProblemIterator) =  decide_convergence(P.M, P.F, P.state, P.options)

"Get last element in the iteration, which is xstar, or throw a warning"
function Base.last(P::ZeroProblemIterator)
#    state.convergence_failed && @warn "The problem failed to converge"
    #    !(state.x_converged || state.f_converged) && @warn "The problem has not converged. Try calling `solve! first"
    P.state.xstar
end

"""
    tracks(P::ZeroProblemIterator)

Show trace of output when `verbose=true` is specified to the problem
"""
function tracks(P::ZeroProblemIterator{M,F,S,O,L}) where {M,F,S,O,L<:Tracks}
     show_trace(P.M, nothing, P.state, P.logger)
end
tracks(P::ZeroProblemIterator) = error("Set verbose=true when specifying the problem to see the tracks")

"""
    solve(fx::ZeroProblem, [p=nothing]; kwargs...)
    solve(fx::ZeroProblem, M, [p=nothing]; kwargs...)
    init(fx::ZeroProblem, [M], [p=nothing];
         verbose=false, tracks=NullTracks(), kwargs...)
    solve!(P::ZeroProblemIterator)

Solve for the zero of a function specified through a  `ZeroProblem` or `ZeroProblemIterator`

The methods involved with this interface are:

* `ZeroProblem`: used to specificy a problem with a function (or functions) and an initial guess
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

The argument `verbose=true` for `init` instructs that steps to be logged; these may be viewed with the method `Roots.tracks` for the iterator.

The iterator interface allows for the creation of hybrid solutions, for example, this is essentially how `Order0` is constructed (`Order0` follows secant steps until a bracket is identified, after which is switches to a bracketing algorithm.)

```
function order0(f, x)
    fx = ZeroProblem(f, x)
    p = init(fx, Roots.Secant())
    xᵢ,st = ϕ = iterate(p)
    while ϕ != nothing
        xᵢ, st = ϕ
        fᵢ₋₁, fᵢ = p.state.fxn0, p.state.fxn1
        if sign(fᵢ₋₁)*sign(fᵢ) < 0 # check for bracket
            x0 = (p.state.xn0, p.state.xn1)
            fx′ = ZeroProblem(f, x0)
            p = init(fx′, Bisection())
            solve!(p)
            break
        end
        ϕ = iterate(p, st)
    end
    return last(p)
end
```

"""
function solve!(P::ZeroProblemIterator)
    log_step(P)            # initial logging
    for _ in P end         # iterate to completion
    last(P)
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
