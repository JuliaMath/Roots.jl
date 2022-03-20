## --------------------------------------------------
## AbstractAcceleratedBisection

# use xatol, xrtol only, but give some breathing room over the strict ones and cap number of steps
function default_tolerances(::AbstractAcceleratedBisection, ::Type{T}, ::Type{S}) where {T,S}
    xatol = eps(real(T))^3 * oneunit(real(T))
    xrtol = 2eps(real(T))  # unitless
    atol = zero(S) * oneunit(real(S))
    rtol = zero(S)
    maxevals = 60
    maxfnevals = typemax(Int)
    strict = false
    (xatol, xrtol, atol, rtol, maxevals, maxfnevals, strict)
end

function init_state(M::AbstractAcceleratedBisection, F, x₀, x₁, fx₀, fx₁)
    (iszero(fx₀) || iszero(fx₁)) && return UnivariateZeroState(x₁, x₀, fx₁, fx₀)
    assert_bracket(fx₀, fx₁)
    a, b, fa, fb = (x₀ < x₁) ? (x₀, x₁, fx₀, fx₁) : (x₁, x₀, fx₁, fx₀)
    UnivariateZeroState(b, a, fb, fa)
end

initial_fncalls(::AbstractAcceleratedBisection) = 2
