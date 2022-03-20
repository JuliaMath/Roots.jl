##################################################
## some means of guarding against large fx when taking a secant step
## TODO: rework this
function steff_step(M::Union{Order5,Order8,Order16}, x::S, fx::T) where {S,T}
    xbar, fxbar = real(x / oneunit(x)), fx / oneunit(fx)
    thresh = max(1, abs(xbar)) * sqrt(eps(one(xbar))) #^(1/2) # max(1, sqrt(abs(x/fx))) * 1e-6

    out = abs(fxbar) <= thresh ? fxbar : sign(fx) * thresh
    x + out * oneunit(x)
end
