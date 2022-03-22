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

incfn(T::AbstractTracks, i=1) = log_fncall(T,i)      # legacy alias
# a tracks object to record tracks
"""
    Tracks(T, S)
    Tracks()


Construct a `Tracks` object used to record the progress of the algorithm.
`T` is the type of function inputs, and `S` is the type of function outputs. They
both default to `Float64`. Note that because this type is not exported, you have to
write `Roots.Tracks()` to construct a `Tracks` object.

By default, no tracking is done while finding a root (with `find_zero`, `solve`, or `fzero`).
To change this, construct a `Tracks` object, and pass it to the keyword argument `tracks`.
This will modify the `Tracks` object, storing the input and function values at each iteration,
along with additional information about the root-finding process.

`Tracker` objects are printed in a nice, easy-to-read
format. Internally either a tuple of `(x,f(x))` pairs or `(aₙ, bₙ)`
pairs are stored, the latter for bracketing methods, as illustrated
in the examples below. (These implementation details may change without notice.)

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
x_0 =  0.0000000000000000,	 fx_0 = -2.0000000000000000
x_1 =  2.0000000000000000,	 fx_1 =  2.0000000000000000
x_2 =  1.0000000000000000,	 fx_2 = -1.0000000000000000
x_3 =  1.3333333333333333,	 fx_3 = -0.2222222222222223
x_4 =  1.4285714285714286,	 fx_4 =  0.0408163265306123
x_5 =  1.4137931034482758,	 fx_5 = -0.0011890606420930
x_6 =  1.4142114384748701,	 fx_6 = -0.0000060072868389
x_7 =  1.4142135626888697,	 fx_7 =  0.0000000008931456
x_8 =  1.4142135623730947,	 fx_8 = -0.0000000000000009

julia> tracker.xfₛ  # stored as (x, f(x)) pairs
9-element Vector{Tuple{Float64, Float64}}:
 (0.0, -2.0)
 (2.0, 2.0)
 (1.0, -1.0)
 (1.3333333333333333, -0.22222222222222232)
 (1.4285714285714286, 0.04081632653061229)
 (1.4137931034482758, -0.0011890606420930094)
 (1.41421143847487, -6.007286838860537e-6)
 (1.4142135626888697, 8.931455575122982e-10)
 (1.4142135623730947, -8.881784197001252e-16)

julia> tracker = Roots.Tracks()
Algorithm has not been run

julia> find_zero(sin, (3, 4), Roots. A42(), tracks=tracker) ≈ π
true

julia> tracker.abₛ # stored as (aₙ, bₙ) pairs
4-element Vector{Tuple{Float64, Float64}}:
 (3.0, 4.0)
 (3.0, 3.5)
 (3.14156188905231, 3.1416247553172956)
 (3.141592653589793, 3.1415926535897936)
```

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
Tracks(T, S) = Tracks(Tuple{T,S}[], Tuple{T,T}[], 0, 0, :not_converged, "",  NaN*zero(T), nothing, nothing)
Tracks(s::AbstractUnivariateZeroState{T,S}) where {T,S} = Tracks(T, S)
Tracks(verbose, tracks, state::AbstractUnivariateZeroState{T,S}) where {T,S} =
    (verbose && isa(tracks, NullTracks)) ? Tracks(T,S) : tracks
Tracks() = Tracks(Float64, Float64) # give default

function log_step(l::Tracks, M::AbstractNonBracketing, state; init=false)
    init && push!(l.xfₛ, (state.xn0, state.fxn0))
    push!(l.xfₛ, (state.xn1, state.fxn1))
    !init  && log_iteration(l, 1)
    nothing
end

log_fncall(l::Tracks, i=1) = (l.fncalls += i; nothing)
log_iteration(l::Tracks, n=1) = (l.steps += n; nothing)
log_message(l::Tracks, msg) = (l.message *= msg; nothing)
log_convergence(l::Tracks, msg) = (l.convergence_flag = msg; nothing)
log_last(l::Tracks, α) = (l.alpha = α; nothing)
log_method(l::Tracks, method) = (l.method = method; nothing)
log_nmethod(l::Tracks, method) = (l.nmethod = method; nothing)

Base.show(io::IO, l::Tracks) = show_trace(io, l.method, l.nmethod, l)


function show_trace(io::IO, method, N, tracks)

    if length(tracks.xfₛ) == 0 && length(tracks.abₛ) == 0
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
                "%s = % 18.16f,\t %s = % 18.16f",
                "x_$(i-1)",
                float(xi),
                "fx_$(i-1)",
                float(fxi)
            )
        )
    end

    # show bracketing
    i₀ = length(s.xfₛ)
    for (i, (a,b)) in enumerate(s.abₛ)
        j = i₀ + i
        println(
            io,
            @sprintf(
                "(%s, %s) = (% 18.16f, % 18.16f )",
                "a_$(j-1)",
                "b_$(j-1)",
                a, b
            )
        )
    end
    println(io, "")

end
