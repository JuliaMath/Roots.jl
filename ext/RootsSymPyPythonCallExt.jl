module RootsSymPyPythonCallExt

using Roots
using SymPyPythonCall

## Allow equations to specify a problem to solve
function Roots.Callable_Function(M::Roots.AbstractUnivariateZeroMethod, f::SymPyPythonCall.Sym, p=nothing)
    if f.is_Equality == true
        f = lhs(f) - rhs(f)
    end
    Roots.Callable_Function(M, lambdify(f), p)
end

function Roots.FnWrapper(f::SymPyPythonCall.Sym)
    if f.is_Equality == true
        f = lhs(f) - rhs(f)
    end
    Roots.FnWrapper(lambdify(f))
end


## allow find_zeros to use symbolic equation
function Roots.find_zeros(f::SymPyPythonCall.Sym, a, b=nothing; kwargs...)
    if f.is_Equality == true
        f = lhs(f) - rhs(f)
    end
    find_zeros(lambdify(f), a, b; kwargs...)
end


end
