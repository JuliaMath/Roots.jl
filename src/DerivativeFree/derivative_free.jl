## Derivative free methods inherit from abstract secant

# init_state(M,F,x) --> call init_state(M,F,x₀,x₁,fx₀, fx₁)
function init_state(M::AbstractSecantMethod, F::Callable_Function, x)
    x₀, x₁ = x₀x₁(x)
    fx₀, fx₁ = first(F(x₀)), first(F(x₁))
    state = init_state(M, F, x₀, x₁, fx₀, fx₁)
end

# initialize from xs, fxs
function init_state(::AbstractSecantMethod, F, x₀, x₁, fx₀, fx₁)
    UnivariateZeroState(x₁, x₀, fx₁, fx₀)
end

initial_fncalls(::AbstractSecantMethod) = 2

##################################################
## Guard against non-robust algorithms
##
## By default, we do this by deciding if we
## should take a secant step instead of the algorithm For example, for
## Steffensen which is quadratically convergent and the Secant method
## which is only superlinear,
## the error, e_{n+1} = x_{n+1} - alpha, may be smaller after a secant
## step than a Steffensen step. (It is only once x_n is close enough
## to alpha that the method is quadratically convergent.
## The Steffensen error is
## Δn+1 = f[x,x+fx, alpha]/f[x, x+fx] * (1 - f[x, alpha]) (x-alpha)^2
##      ≈ f''/(2f') * ( 1 + f') Δn^2
## The Secant error is
## Δn+1 = f[x,x_{-1},alpha] / f[x,x_{-1}] * (x-alpha) * (x_{-1} - alpha)
##      ≈  f''/(2f')  Δn ⋅ Δn-1
## The ratio is ≈ (1 + f')(Δn / Δn-1)
## It seems reasonable, that a Steffensen step is preferred when
## the ratio satisfies -1 < (1+f') ⋅ Δn /Δn-1 < 1
## We could use f' ~ fp = (fx1-fx0)/(x1-x0); but our proxy for
## Δn/Δn-1 is problematic, as we don't know alpha, and using xn-x_{n-1}
## can be an issue when only x1 and not x0 is specified. This needs
## working around.
##
## Instead, as Steffensen is related to Newton as much as
## (f(x+fx) - fx)/fx  ≈ f'(x), we take a Steffensen step if |fx|
## is small enough. For this we use |fx| <= x/1000; which
## seems to work reasonably well over several different test cases.

@inline function do_guarded_step(
    M::AbstractSecantMethod,
    o::AbstractUnivariateZeroState{T,S},
) where {T,S}
    x, fx = o.xn1, o.fxn1
    1000 * abs(fx) > max(oneunit(S), abs(x) * oneunit(S) / oneunit(T)) * one(T)
end

# check if we should guard against step for method M; call N if yes, P if not
function update_state_guarded(
    M::AbstractSecantMethod,
    N::AbstractUnivariateZeroMethod,
    P::AbstractUnivariateZeroMethod,
    fs,
    o,
    options,
    l=NullTracks(),
)
    if do_guarded_step(M, o)
        return update_state(N, fs, o, options, l)
    else
        update_state(P, fs, o, options, l)
    end
end
