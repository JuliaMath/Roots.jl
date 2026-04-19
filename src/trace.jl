### --------------------------------------------------

## Tracks (for logging actual steps)
## when no logging this should get optimized out
## when logging, this allocates

abstract type AbstractTracks end

"""
    Roots.Tracks()

A `Tracks` instance is used to record the progress of an algorithm.
Note that because this type is not exported, you have to
call `Roots.Tracks()` to construct a `Tracks` object.

By default, no tracking is done while finding a root.
To change this, construct a `Tracks` object, and pass it to the keyword argument `tracks`.
This will modify the `Tracks` object, storing the input and function values at each iteration,
along with additional information about the root-finding process.

`Tracks` objects are shown in an easy-to-read
format. Internally either a tuple of `(x,f(x))` pairs or `(aₙ, bₙ)`
pairs are stored, the latter for bracketing methods. (These
implementation details may change without notice.).


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
x₁ =  0                   	 fx₁ = -2
x₂ =  2                   	 fx₂ =  2
x₃ =  1                   	 fx₃ = -1
x₄ =  1.3333333333333333  	 fx₄ = -0.22222222222222232
x₅ =  1.4285714285714286  	 fx₅ =  0.04081632653061229
x₆ =  1.4137931034482758  	 fx₆ = -0.0011890606420930094
x₇ =  1.4142114384748701  	 fx₇ = -6.0072868388605372e-06
x₈ =  1.4142135626888697  	 fx₈ =  8.9314555751229818e-10
x₉ =  1.4142135623730947  	 fx₉ = -8.8817841970012523e-16

julia> tracker = Roots.Tracks()
Algorithm has not been run

julia> find_zero(sin, (3, 4), Roots.A42(), tracks=tracker) ≈ π
true

julia> tracker
Results of univariate zero finding:

* Converged to: 3.141592653589793
* Algorithm: A42()
* iterations: 3
* function evaluations ≈ 9
* stopped as x_n ≈ x_{n-1} using atol=xatol, rtol=xrtol

Trace:
(a₁, b₁) = (  3                  ,  4                   )
(a₂, b₂) = (  3                  ,  3.157162792479947   )
(a₃, b₃) = (  3.1415926144917452 ,  3.1415926926910007  )
(a₄, b₄) = (  3.1415926535897931 ,  3.141592653589794   )
```

!!! note
    As designed, the `Tracks` object may not record actions taken
    while the state object is initialized. An example is the default
    bisection algorithm where an initial midpoint is found to ensure
    the bracket does not straddle ``0``.

"""
mutable struct Tracks <: AbstractTracks
    h
    steps::Int
    fncalls::Int
    convergence_flag::Symbol
    message::String
    alpha
    method
    nmethod
end

# mimic MVHistory from ValueHistories
# we only need a fraction of the feature set and avoid dependencies
struct MVHistory
    d::Dict{Symbol,Any}
end

function MVHistory()
    d = Dict{Symbol,Any}()
    MVHistory(d)
end

Base.haskey(h::MVHistory, k::Symbol) = haskey(h.d, k)
function Base.push!(h::MVHistory, k::Symbol, i, val)
    if !haskey(h, k)
        h.d[k] = (Any[i], Any[val])
    else
        inds, vals = h.d[k]
        push!(inds, i)
        push!(vals, val)
    end
end
Base.get(h::MVHistory, k::Symbol) = h.d[k]
Base.length(h::MVHistory, k::Symbol) = length(first(h.d[k]))

Tracks() = Tracks(MVHistory(), 0, 0, :algorithm_not_run, "", NaN, nothing, nothing)

# default for no logging
struct NullTracks <: AbstractTracks end

# API

# how many function evaluations in init_state
# this is a worst case estimate leading to the function call count being an upper bound only
initial_fncalls(M::AbstractUnivariateZeroState) = @warn "initial_fncalls fix $M"

# For NullTracks do nothing
log_step(s::NullTracks, M, x; init=false) = nothing  # log a step (x,f(x)) or (a,b)
log_iteration(::NullTracks, n=1) = nothing           # log an iteration (call to update state)
log_fncall(::NullTracks, n=1) = nothing              # add a function call (aliased to incfn); only one called in update_step
log_message(::NullTracks, msg) = nothing             # append to a message
log_convergence(::NullTracks, msg) = nothing         # flag for convergence
log_last(::NullTracks, state) = nothing              # log α
log_method(::NullTracks, method) = nothing           # record M
log_nmethod(::NullTracks, method) = nothing          # record N (if hybrid)

incfn(T::AbstractTracks, i=1) = log_fncall(T, i)      # legacy alias

# for Tracks
## There are specializations on M in bracketing, AlefeldPotraShi, Lith, and Newton,
function log_step(l::Tracks, M::AbstractNonBracketingMethod, x; init=false)
    h, 𝑀 = l.h, nameof(typeof(M))
    init && push!(h, 𝑀, 1, (x.xn0, x.fxn0))
    n = haskey(h, 𝑀) ? length(h, 𝑀) : 1
    push!(h, 𝑀, n + 1, (x.xn1, x.fxn1))
    !init && log_iteration(l, 1)
    nothing
end

log_fncall(l::Tracks, i=1) = (l.fncalls += i; nothing)
log_iteration(l::Tracks, n=1) = (l.steps += n; nothing)
function log_message(l::Tracks, msg)
    cur = l.message
    l.message = join(cur, msg)
    nothing
end

log_convergence(l::Tracks, msg) = (l.convergence_flag=msg; nothing)
log_last(l::Tracks, α) = (l.alpha=α; nothing)
log_method(l::Tracks, method) = (l.method=method; nothing)
log_nmethod(l::Tracks, method) = (l.nmethod=method; nothing)

# reset tracker
Base.empty!(l::NullTracks) = nothing
function Base.empty!(l::Tracks)
    l.h = MVHistory()
    l.steps = 0
    l.fncalls = 0
    l.convergence_flag = :algorithm_not_run
    l.message = ""
    l.alpha = NaN
    l.method = l.nmethod = nothing
    nothing
end

## --- display tracks
Base.show(io::IO, l::NullTracks) = nothing

Base.show(io::IO, l::Tracks) = show_trace(io, l.method, l.nmethod, l)

function show_trace(io::IO, M, N, tracks)
    flag = tracks.convergence_flag
    if flag == :algorithm_not_run
        print(io, "Algorithm has not been run")
        return nothing
    end

    𝑀, 𝑁 = nameof(typeof(M)), nameof(typeof(N))
    n = haskey(tracks.h, 𝑁) ? length(tracks.h, 𝑀) : 0
    converged = !isnan(tracks.alpha)
    println(io, "Results of univariate zero finding:\n")
    if converged
        println(io, "* Converged to: $(tracks.alpha)")
        if N === nothing || n == 0
            println(io, "* Algorithm: $(M)")
        else
            println(io, "* Algorithm: $(M); finished with bracketing method $N")
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
        println(io, "* Algorithm $(M)")
    end
    println(io, "")
    println(io, "Trace:")
    show_tracks(io, tracks, M)
    !isnothing(N) && show_tracks(io, tracks, N)
end

# XXX of @printf here triggers a possible error flag with JET.
function show_tracks(io::IO, s::Tracks, M::AbstractUnivariateZeroMethod)
    # show (x,f(x))
    𝑀 = nameof(typeof(M))
    ind, xf = get(s.h, 𝑀)
    for (i, (xi, fxi)) in zip(ind, xf)
        𝑖 = sprint(io -> unicode_subscript(io, i))
        if !(isa(xi, Complex) || isa(fxi, Complex))
            @printf(
                io,
                "%s%s = % -20.17g \t %s%s = % -20.17g\n",
                "x",
                𝑖,
                float(xi),
                "fx",
                𝑖,
                float(fxi)
            )
        else
            @printf(
                io,
                "%s%s = (% -20.17g, % -20.17g),\t %s%s = (% -20.17g, % -20.17g)\n",
                "x",
                𝑖,
                real(xi),
                imag(xi),
                "fx",
                𝑖,
                real(fxi),
                imag(fxi)
            )
        end
    end
end

function show_tracks(io::IO, s::Tracks, M::AbstractBracketingMethod)
    # show (a,b)
    h = s.h
    𝑀 = nameof(typeof(M))
    𝑁 = nameof(typeof(s.nmethod))
    i₀ = haskey(h.d, 𝑁) ? length(h, 𝑀) : 0
    ind, ab = get(s.h, 𝑀)

    for (i, (a, b)) in zip(ind, ab)
        j = i₀ + i
        𝑗 = sprint(io -> unicode_subscript(io, j))
        @printf(io, "(%s%s, %s%s) = ( % -20.17g, % -20.17g )\n", "a", 𝑗, "b", 𝑗, a, b)
    end
end

## --- these could be deleted as methods ...
#=
find_zerov(f, x, M; kwargs...)

Run `find_zero` return a `Tracks` object, not the value, which can be extracted via the `last` method.
=#
function find_zerov(f, x, M; tracks=nothing, kwargs...)
    tracks′ = isnothing(tracks) ? Roots.Tracks() : tracks
    Z = init(ZeroProblem(f, x), M; tracks=tracks′, kwargs...)
    solve!(Z)
    tracks′
end
find_zerov(f, x; tracks=nothing, kwargs...) =
    find_zerov(f, x, find_zero_default_method(x); tracks=tracks, kwargs...)
