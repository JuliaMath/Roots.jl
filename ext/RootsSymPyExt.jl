module RootsSymPyExt

using Roots
using SymPy

## Allow equations to specify a problem to solve
function Roots.Callable_Function(M::Roots.AbstractUnivariateZeroMethod, f::SymPy.Sym, p=nothing)
    if f.is_Equality == true
        f = lhs(f) - rhs(f)
    end
    Roots.Callable_Function(M, lambdify(f), p)
end


end
