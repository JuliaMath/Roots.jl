### --------------------------------------------------

## Tracks (for logging actual steps)
## when no logging this should get optimized out
## when logging, this allocates

abstract type AbstractTracks end
struct NullTracks <: AbstractTracks end

# logging api

# how many function evaluations in init_state
# this is a worst case estimate leading to the function call count being an upper bound only
initial_fncalls(M::AbstractUnivariateZeroState) = @warn "initial_fncalls fix $M"

log_step(s::NullTracks, M, x; init=false) = nothing  # log a step (x,f(x)) or (a,b)
log_iteration(::NullTracks, n=1) = nothing           # log an iteration (call to update state)
log_fncall(::NullTracks, n=1) = nothing              # add a function call (aliased to incfn); only one called in update_step
log_message(::NullTracks, msg) = nothing             # append to a message
log_convergence(::NullTracks, msg) = nothing         # flag for convergence
log_last(::NullTracks, state) = nothing              # log α
log_method(::NullTracks, method) = nothing           # record M
log_nmethod(::NullTracks, method) = nothing          # record N (if hybrid)

incfn(T::AbstractTracks, i=1) = log_fncall(T, i)      # legacy alias
# a tracks object to record tracks
"""
    Tracks(T, S)
    Tracks()


Construct a `Tracks` object used to record the progress of the algorithm.
`T` is the type of function inputs, and `S` is the type of function outputs. They
both default to `Float64`. Note that because this type is not exported, you have to
write `Roots.Tracks()` to construct a `Tracks` object.

By default, no tracking is done while finding a root.
To change this, construct a `Tracks` object, and pass it to the keyword argument `tracks`.
This will modify the `Tracks` object, storing the input and function values at each iteration,
along with additional information about the root-finding process.

`Tracks` objects are shown in an easy-to-read
format. Internally either a tuple of `(x,f(x))` pairs or `(aₙ, bₙ)`
pairs are stored, the latter for bracketing methods. (These
implementation details may change without notice.) The methods
`empty!`, to reset the `Tracks` object; `get`, to get the tracks;
`last`, to get the value converted to, may be of interest.

If you only want to print the information, but you don't need it later, this can conveniently be
done by passing `verbose=true` to the root-finding function. This will not
effect the return value, which will still be the root of the function.


## Examples

```jldoctest Tracks
julia> using Roots

julia> f(x) = x^2-2
f (generic function with 1 method)

julia> tracker = Roots.Tracks()
Algorithm has not been run

julia> find_zero(f, (0, 2), Roots.Secant(), tracks=tracker) ≈ √2
true

julia> tracker
Results of univariate zero finding:

* Converged to: 1.4142135623730947
* Algorithm: Secant()
* iterations: 7
* function evaluations ≈ 9
* stopped as |f(x_n)| ≤ max(δ, |x|⋅ϵ) using δ = atol, ϵ = rtol

Trace:
x₁ = 0,	 fx₁ = -2
x₂ = 2,	 fx₂ = 2
x₃ = 1,	 fx₃ = -1
x₄ = 1.3333333333333333,	 fx₄ = -0.22222222222222232
x₅ = 1.4285714285714286,	 fx₅ = 0.04081632653061229
x₆ = 1.4137931034482758,	 fx₆ = -0.0011890606420930094
x₇ = 1.4142114384748701,	 fx₇ = -6.0072868388605372e-06
x₈ = 1.4142135626888697,	 fx₈ = 8.9314555751229818e-10
x₉ = 1.4142135623730947,	 fx₉ = -8.8817841970012523e-16

julia> empty!(tracker)  # resets

julia> find_zero(sin, (3, 4), Roots.A42(), tracks=tracker) ≈ π
true

julia> get(tracker)
4-element Vector{NamedTuple{names, Tuple{Float64, Float64}} where names}:
 (a = 3.0, b = 4.0)
 (a = 3.0, b = 3.157162792479947)
 (a = 3.141592614491745, b = 3.1415926926910007)
 (a = 3.141592653589793, b = 3.141592653589794)

julia> last(tracker)
3.141592653589793
```

!!! note
    As designed, the `Tracks` object may not record actions taken
    while the state object is initialized. An example is the default
    bisection algorithm where an initial midpoint is found to ensure
    the bracket does not straddle ``0``.

"""
mutable struct Tracks{T,S} <: AbstractTracks
    xfₛ::Vector{Tuple{T,S}} # (x,f(x))
    abₛ::Vector{Tuple{T,T}} # (aᵢ, bᵢ)
    steps::Int
    fncalls::Int
    convergence_flag::Symbol
    message::String
    alpha::T
    method
    nmethod
end
Tracks(T, S) = Tracks(
    Tuple{T,S}[],
    Tuple{T,T}[],
    0,
    0,
    :algorithm_not_run,
    "",
    nan(T),
    nothing,
    nothing,
)
Tracks(s::AbstractUnivariateZeroState{T,S}) where {T,S} = Tracks(T, S)
Tracks(verbose, tracks, state::AbstractUnivariateZeroState{T,S}) where {T,S} =
    (verbose && isa(tracks, NullTracks)) ? Tracks(T, S) : tracks
Tracks() = Tracks(Float64, Float64) # give default

function log_step(l::Tracks, M::AbstractNonBracketingMethod, state; init=false)
    init && push!(l.xfₛ, (state.xn0, state.fxn0))
    push!(l.xfₛ, (state.xn1, state.fxn1))
    !init && log_iteration(l, 1)
    nothing
end

log_fncall(l::Tracks, i=1) = (l.fncalls += i; nothing)
log_iteration(l::Tracks, n=1) = (l.steps += n; nothing)
log_message(l::Tracks, msg) = (l.message *= msg; nothing)
log_convergence(l::Tracks, msg) = (l.convergence_flag = msg; nothing)
log_last(l::Tracks, α) = (l.alpha = α; nothing)
log_method(l::Tracks, method) = (l.method = method; nothing)
log_nmethod(l::Tracks, method) = (l.nmethod = method; nothing)

# copy some API from ValueHistories
Base.first(l::AbstractTracks) = (@warn "No tracking information was kept"; nothing)
function Base.first(l::Tracks)
    l.convergence_flag == :algorithm_not_run && error("Algorithm not run")
    !isempty(l.xfₛ) && return (x₁=l.xfₛ[1][1], f₁=l.xfₛ[1][2])
    return (a₁=l.abₛ[1][1], b₁=l.abₛ[1][2])
end
# return last value of algorithm
Base.last(l::AbstractTracks) = (@warn "No tracking information was kept"; nothing)
function Base.last(l::Tracks)
    convergence_flag = l.convergence_flag
    α = l.alpha
    if convergence_flag == :algorithm_not_run
        @warn "The algorithm has not run"
    end
    α
end

# Returns all available observations.
Base.get(l::NullTracks) = (@warn "No tracking information was kept"; nothing)
function Base.get(l::Tracks)
    xf = [(xn=xᵢ, fn=fᵢ) for (xᵢ, fᵢ) in l.xfₛ]
    ab = [(a=min(u, v), b=max(u, v)) for (u, v) in l.abₛ]
    vcat(xf, ab)
end

# reset tracker
Base.empty!(l::NullTracks) = nothing
function Base.empty!(l::Tracks{T,S}) where {T,S}
    empty!(l.xfₛ)
    empty!(l.abₛ)
    l.steps = 0
    l.fncalls = 0
    l.convergence_flag = :algorithm_not_run
    l.message = ""
    l.alpha = nan(T)
    l.method = l.nmethod = nothing
    nothing
end

Base.show(io::IO, l::Tracks) = show_trace(io, l.method, l.nmethod, l)

function show_trace(io::IO, method, N, tracks)
    flag = tracks.convergence_flag
    if flag == :algorithm_not_run
        print(io, "Algorithm has not been run")
        return nothing
    end

    converged = !isnan(tracks.alpha)
    println(io, "Results of univariate zero finding:\n")
    if converged
        println(io, "* Converged to: $(tracks.alpha)")
        if N === nothing || length(tracks.abₛ) == 0
            println(io, "* Algorithm: $(method)")
        else
            println(io, "* Algorithm: $(method); finished with bracketing method $N")
        end
        println(io, "* iterations: $(tracks.steps)")
        println(io, "* function evaluations ≈ $(tracks.fncalls)")
        tracks.convergence_flag == :x_converged &&
            println(io, "* stopped as x_n ≈ x_{n-1} using atol=xatol, rtol=xrtol")
        tracks.convergence_flag == :f_converged &&
            tracks.message == "" &&
            println(io, "* stopped as |f(x_n)| ≤ max(δ, |x|⋅ϵ) using δ = atol, ϵ = rtol")
        tracks.convergence_flag == :exact_zero && println(io, "* stopped as f(x_n) = 0")
        tracks.message != "" && println(io, "* Note: $(tracks.message)")
    else
        println(io, "* Convergence failed: $(tracks.message)")
        println(io, "* Algorithm $(method)")
    end
    println(io, "")
    println(io, "Trace:")
    show_tracks(io, tracks, method)
end

function show_tracks(io::IO, s::Tracks, M::AbstractUnivariateZeroMethod)

    # show (x,f(x))
    for (i, (xi, fxi)) in enumerate(s.xfₛ)
        println(
            io,
            @sprintf(
                "%s%s = %.17g,\t %s%s = %.17g",
                "x",
                sprint(io -> unicode_subscript(io, i)),
                float(xi),
                "fx",
                sprint(io -> unicode_subscript(io, i)),
                float(fxi)
            )
        )
    end

    # show bracketing
    i₀ = length(s.xfₛ)
    for (i, (a, b)) in enumerate(s.abₛ)
        j = i₀ + i
        println(
            io,
            @sprintf(
                "(%s%s, %s%s) = ( %.17g, %.17g )",
                "a",
                sprint(io -> unicode_subscript(io, j - 1)),
                "b",
                sprint(io -> unicode_subscript(io, j - 1)),
                a,
                b
            )
        )
    end
    println(io, "")
end

# support for complex values
# Issue 336. (Could DRY this up...)
function show_tracks(
    io::IO,
    s::Roots.Tracks{T,S},
    M::Roots.AbstractUnivariateZeroMethod,
) where {T<:Complex,S<:Complex}

    # show (x,f(x))
    for (i, (xi, fxi)) in enumerate(s.xfₛ)
        println(
            io,
            @sprintf(
                "%s%s = (%.17g, %.17g),\t %s%s = (%.17g, %.17g)",
                "x",
                sprint(io -> Roots.unicode_subscript(io, i)),
                real(xi),
                imag(xi),
                "fx",
                sprint(io -> Roots.unicode_subscript(io, i)),
                real(fxi),
                imag(fxi)
            )
        )
    end

    # show bracketing
    i₀ = length(s.xfₛ)
    for (i, (a, b)) in enumerate(s.abₛ)
        j = i₀ + i
        println(
            io,
            @sprintf(
                "(%s%s, %s%s) = ( %.17g, %.17g )",
                "a",
                sprint(io -> unicode_subscript(io, j - 1)),
                "b",
                sprint(io -> unicode_subscript(io, j - 1)),
                a,
                b
            )
        )
    end
    println(io, "")
end

## needs better name, but is this useful?
"""
    find_zerov(f, x, M; kwargs...)

Run `find_zero` return a `Tracks` object, not the value, which can be extracted via the `last` method.
"""
function find_zerov(f, x, M; verbose=nothing, kwargs...)
    Z = init(ZeroProblem(f, x), M; verbose=true, kwargs...)
    solve!(Z)
    Z.logger
end
find_zerov(f, x; verbose=nothing, kwargs...) =
    find_zerov(f, x, find_zero_default_method(x); kwargs...)
