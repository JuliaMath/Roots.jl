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

# return f(x); (f(x), f(x)/f'(x)); *or* f(x), (f(x)/f'(x), f'(x)/f''(x), ...) # so N=1, 2 are special cased
# Callable_Function(output_arity, input_arity, F, p)
# First handle: x -> (f,f/f', f'/f'', ...)
(F::Callable_Function{Val{1},Val{false},𝑭,Nothing})(x) where {𝑭} = first(F.f(x))
(F::Callable_Function{Val{1},Val{false},𝑭,P})(x) where {𝑭,P} = first(F.f(x, F.p))

(F::Callable_Function{Val{2},Val{false},𝑭,Nothing})(x) where {𝑭} = F.f(x)[1:2]
(F::Callable_Function{Val{2},Val{false},𝑭,P})(x) where {𝑭,P} = F.f(x, F.p)[1:2]

# N ≥ 3 returns (f, (...))
function (F::Callable_Function{Val{N},Val{false},𝑭,Nothing})(x) where {N,𝑭}
    fs = F.f(x)
    fs[1], ntuple(i -> fs[i + 1], Val(N - 1))
end
function (F::Callable_Function{Val{N},Val{false},𝑭,P})(x) where {N,𝑭,P}
    fs = F.f(x, F.p)
    fs[1], ntuple(i -> fs[i + 1], Val(N - 1))
end

## f is specified as a tuple (f,f',f'', ...)
## N =1  return f(x)
(F::Callable_Function{Val{1},Val{true},𝑭,Nothing})(x) where {𝑭} = first(F.f)(x)
(F::Callable_Function{Val{1},Val{true},𝑭,P})(x) where {𝑭,P} = first(F.f)(x, F.p)
## N=2 return (f(x), f(x)/f'(x))
function (F::Callable_Function{Val{2},Val{true},𝑭,Nothing})(x) where {𝑭}
    f, f′ = (F.f[1])(x), (F.f[2])(x)
    f, f / f′
end
function (F::Callable_Function{Val{2},Val{true},𝑭,P})(x) where {𝑭,P}
    f, f′ = (F.f[1])(x, F.p), (F.f[2])(x, F.p)
    f, f / f′
end

## For N ≥ 3 we return (f, (f/f', f'/f'', ...);
## Pay no attention to this code; we hand write a bunch, as the
## general formula later runs more slowly.
function (F::Callable_Function{Val{3},Val{true},𝑭,Nothing})(x) where {𝑭}
    f, f′, f′′ = (F.f[1])(x), (F.f[2])(x), (F.f[3])(x)
    f, (f / f′, f′ / f′′)
end
function (F::Callable_Function{Val{3},Val{true},𝑭,P})(x) where {𝑭,P}
    f, f′, f′′ = (F.f[1])(x, F.p), (F.f[2])(x, F.p), (F.f[3])(x, F.p)
    f, (f / f′, f′ / f′′)
end

function (F::Callable_Function{Val{4},Val{true},𝑭,Nothing})(x) where {𝑭}
    f, f′, f′′, f′′′ = (F.f[1])(x), (F.f[2])(x), (F.f[3])(x), (F.f[4])(x)
    f, (f / f′, f′ / f′′, f′′ / f′′′)
end
function (F::Callable_Function{Val{4},Val{true},𝑭,P})(x) where {𝑭,P}
    f, f′, f′′, f′′′ =
        (F.f[1])(x, F.p), (F.f[2])(x, F.p), (F.f[3])(x, F.p), (F.f[4])(x, F.p)
    𝐓 = eltype(f / f′)
    f, NTuple{3,𝐓}((f / f′, f′ / f′′, f′′ / f′′′))
end

function (F::Callable_Function{Val{5},Val{true},𝑭,Nothing})(x) where {𝑭}
    f, f′, f′′, f′′′, f′′′′ =
        (F.f[1])(x), (F.f[2])(x), (F.f[3])(x), (F.f[4])(x), (F.f[5])(x)
    f, (f / f′, f′ / f′′, f′′ / f′′′, f′′′ / f′′′′)
end
function (F::Callable_Function{Val{5},Val{true},𝑭,P})(x) where {𝑭,P}
    f, f′, f′′, f′′′, f′′′′ = (F.f[1])(x, F.p),
    (F.f[2])(x, F.p),
    (F.f[3])(x, F.p),
    (F.f[4])(x, F.p),
    (F.f[5])(x, F.p)
    f, (f / f′, f′ / f′′, f′′ / f′′′, f′′′ / f′′′′)
end

function (F::Callable_Function{Val{6},Val{true},𝑭,Nothing})(x) where {𝑭}
    f, f′, f′′, f′′′, f′′′′, f′′′′′ =
        (F.f[1])(x), (F.f[2])(x), (F.f[3])(x), (F.f[4])(x), (F.f[5])(x), (F.f[6])(x)
    f, (f / f′, f′ / f′′, f′′ / f′′′, f′′′ / f′′′′, f′′′′ / f′′′′′)
end
function (F::Callable_Function{Val{6},Val{true},𝑭,P})(x) where {𝑭,P}
    f, f′, f′′, f′′′, f′′′′, f′′′′′ = (F.f[1])(x, F.p),
    (F.f[2])(x, F.p),
    (F.f[3])(x, F.p),
    (F.f[4])(x, F.p),
    (F.f[5])(x, F.p),
    (F.f[6])(x, F.p)
    f, (f / f′, f′ / f′′, f′′ / f′′′, f′′′ / f′′′′, f′′′′ / f′′′′′)
end

# faster with the above written out, should generate them...
function (F::Callable_Function{Val{𝐍},Val{true},𝑭,Nothing})(x) where {𝐍,𝑭}
    fs = ntuple(i -> F.f[i](x), Val(𝐍))
    first(fs), ntuple(i -> fs[i] / fs[i + 1], Val(𝐍 - 1))
end

function (F::Callable_Function{Val{𝐍},Val{true},𝑭,P})(x) where {𝐍,𝑭,P}
    fs = ntuple(i -> F.f[i](x, F.p), Val(𝐍))
    first(fs), ntuple(i -> fs[i] / fs[i + 1], Val(𝐍 - 1))
end
