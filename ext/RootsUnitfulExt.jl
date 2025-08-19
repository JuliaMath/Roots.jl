module RootsUnitfulExt
using Unitful
using Roots

function Roots.find_zero_default_method(x0::Tuple{<:Quantity,<:Quantity})
    T = eltype(ustrip.(float.(Roots._extrema(x0))))
    T <: Union{Float16,Float32,Float64} ? Roots.Bisection() : Roots.A42()
end

end # module
