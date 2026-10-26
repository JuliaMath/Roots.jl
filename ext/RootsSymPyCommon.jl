# shared by RootsSymPyExt and RootsSymPyPythonCallExt
# expects `Sym`, `lhs`, `rhs`, and `lambdify` to be defined by the including module

# an equation `lhs ~ rhs` is solved as `lhs - rhs = 0`
_as_expr(f::Sym) = f.is_Equality == true ? lhs(f) - rhs(f) : f

## Allow equations to specify a problem to solve
function Roots.Callable_Function(M::Roots.AbstractUnivariateZeroMethod, f::Sym, p=nothing)
    Roots.Callable_Function(M, lambdify(_as_expr(f)), p)
end

Roots.FnWrapper(f::Sym) = Roots.FnWrapper(lambdify(_as_expr(f)))

## allow find_zeros to use symbolic equation
function Roots.find_zeros(f::Sym, a, b=nothing; kwargs...)
    find_zeros(lambdify(_as_expr(f)), a, b; kwargs...)
end
