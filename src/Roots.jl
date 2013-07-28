module Roots
using Polynomial

export fzero, newton, halley
export thukral
export multroot
export D, D2

include("fzero.jl")
include("adiff.jl")
include("newton.jl")
include("thukral.jl")
include("multroot.jl")

end

