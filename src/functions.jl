## Wrapper for functions to abstract out f, (f,f'), ...

# indicate if we expect f() to return one or multiple values (e.g. Newton)
fn_argout(::AbstractUnivariateZeroMethod) = 1

# A hacky means to call a function so that parameters can be passed as desired
# and the correct number of outputs are computed
struct Callable_Function{Single,Tup,F,P}
    f::F
    p::P
end

function Callable_Function(M::AbstractUnivariateZeroMethod, f, p=nothing)
    Single = Val{fn_argout(M)}
    Tup = Val{isa(f, Tuple)}
    F = typeof(f)
    p′ = ismissing(p) ? nothing : p
    P′ = typeof(p′)
    Callable_Function{Single,Tup,F,P′}(f, p′)
end

function Callable_Function(M::AbstractUnivariateZeroMethod, F::Callable_Function, p=F.p)
    Callable_Function(M, F.f, p)
end

# call g(x) or g(x, p), depending on whether F has a parameter
_apply(::Callable_Function{S,T,𝑭,Nothing}, g, x) where {S,T,𝑭} = g(x)
_apply(F::Callable_Function, g, x) = g(x, F.p)

# return f(x); (f(x), f(x)/f'(x)); *or* f(x), (f(x)/f'(x), f'(x)/f''(x), ...) # so N=1, 2 are special cased
# Callable_Function(output_arity, input_arity, F, p)
# First handle: x -> (f,f/f', f'/f'', ...)
(F::Callable_Function{Val{1},Val{false}})(x) = first(_apply(F, F.f, x))

(F::Callable_Function{Val{2},Val{false}})(x) = _apply(F, F.f, x)[1:2]

# N ≥ 3 returns (f, (...))
function (F::Callable_Function{Val{N},Val{false}})(x) where {N}
    fs = _apply(F, F.f, x)
    fs[1], ntuple(i -> fs[i + 1], Val(N - 1))
end

## f is specified as a tuple (f,f',f'', ...)
## N =1  return f(x)
(F::Callable_Function{Val{1},Val{true}})(x) = _apply(F, first(F.f), x)
## N=2 return (f(x), f(x)/f'(x))
function (F::Callable_Function{Val{2},Val{true}})(x)
    f, f′ = _apply(F, F.f[1], x), _apply(F, F.f[2], x)
    f, f / f′
end

## For N ≥ 3 we return (f, (f/f', f'/f'', ...);
## Pay no attention to this code; we hand write a bunch, as the
## general formula later runs more slowly.
function (F::Callable_Function{Val{3},Val{true}})(x)
    f, f′, f′′ = _apply(F, F.f[1], x), _apply(F, F.f[2], x), _apply(F, F.f[3], x)
    f, (f / f′, f′ / f′′)
end

function (F::Callable_Function{Val{4},Val{true}})(x)
    f, f′, f′′, f′′′ = _apply(F, F.f[1], x),
    _apply(F, F.f[2], x),
    _apply(F, F.f[3], x),
    _apply(F, F.f[4], x)
    f, (f / f′, f′ / f′′, f′′ / f′′′)
end

function (F::Callable_Function{Val{5},Val{true}})(x)
    f, f′, f′′, f′′′, f′′′′ = _apply(F, F.f[1], x),
    _apply(F, F.f[2], x),
    _apply(F, F.f[3], x),
    _apply(F, F.f[4], x),
    _apply(F, F.f[5], x)
    f, (f / f′, f′ / f′′, f′′ / f′′′, f′′′ / f′′′′)
end

function (F::Callable_Function{Val{6},Val{true}})(x)
    f, f′, f′′, f′′′, f′′′′, f′′′′′ = _apply(F, F.f[1], x),
    _apply(F, F.f[2], x),
    _apply(F, F.f[3], x),
    _apply(F, F.f[4], x),
    _apply(F, F.f[5], x),
    _apply(F, F.f[6], x)
    f, (f / f′, f′ / f′′, f′′ / f′′′, f′′′ / f′′′′, f′′′′ / f′′′′′)
end

# faster with the above written out, should generate them...
function (F::Callable_Function{Val{𝐍},Val{true}})(x) where {𝐍}
    fs = ntuple(i -> _apply(F, F.f[i], x), Val(𝐍))
    first(fs), ntuple(i -> fs[i] / fs[i + 1], Val(𝐍 - 1))
end
