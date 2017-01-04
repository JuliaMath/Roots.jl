## Polynomial root finder for polynomials with multiple roots
##
## Based on "Computing multiple roots of inexact polynomials"
## http://www.neiu.edu/~zzeng/mathcomp/zroot.pdf
## Author: Zhonggang Zeng 
## Journal: Math. Comp. 74 (2005), 869-903 
##
## Zeng has a MATLAB package `multroot`, from which this name is derived.
## Basic idea is
## 1) for polynomial p we do gcd decomposition p = u * v; p' = u * w. Then roots(v) are the roots without multiplicities.
## 2) can repeat with u to get multiplicities.
##
## This is from Gauss, as explained in paper. Zeng shows how to get u,v,w when the polynomials
## are inexact due to floating point approximations or even model error. This is done in his 
## algorithm II.
## 3) Zeng's algorithm I (pejroot) uses the pejorative manifold of Kahan and Gauss-Newton to
## improve the root estimates from algorithm II (roots(v)). The pejorative manifold is defined by
## the multiplicities l and is operationalized in evalG and evalJ from Zeng's paper.   

using Polynomials
import Polynomials: degree


## map monic(p) to a point in C^n
## p = 1x^n + a1x^n-1 + ... + an_1 x + an -> (a1,a2,...,an)
function p2a(p::Poly) 
    p = monic(p)
    rcoeffs(p)[2:end]
end

## get value of gl(z). From p16
function evalG(zs::Vector, ls::Vector)
    length(zs) == length(ls) || throw("Length mismatch")

    s = prod([poly([z])^l for (z,l) in zip(zs, ls)])  # \prod (x-z_i)^l_i
    p2a(s)
#    rcoeffs(s)[2:end]
end

## get jacobian J_l(z), p16
function evalJ(zs::Vector, ls::Vector)
    length(zs) == length(ls) || throw("Length mismatch")
    m = length(zs)

    u = prod([poly([z])^(l-1) for (z,l) in zip(zs, ls)]) ## Pi (1-z)^(l-1)

    J = zeros(eltype(zs), sum(ls), m)
    for j in 1:m
        s = -ls[j] * u
        for i in 1:m
            if i != j
                s = s * poly([zs[i]])
            end
        end
        J[:,j] = rcoeffs(s)
    end
    J
end

## Gauss-Newton iteration to solve weighted least squares problem
## G_l(z) = a, where a is related to monic version of polynomial p
## l is known multiplicity structure of polynomial p = (x-z1)^l1 * (x-z2)^l2 * ... * (x-zn)^ln
## Algorithm I, p17
function pejroot(p::Poly, z0::Vector, l::Vector{Int};
                 wts::(@compat Union{Vector, Void})=nothing, # weight vector
                 tol = 1e-8,
                 maxsteps = 100
                      )
    
    a = p2a(p) #rcoeffs(monic(p))[2:end] # an_1, an_2, ..., a2, a1, a0

    if wts == nothing
        @compat wts = map(u -> min(1, 1/abs.(u)), a)
    end
    W = diagm(wts)

    ## Solve WJ Δz = W(Gl(z) - a) in algorithm I
    G(z) = (evalG(z, l) - a)
    update(z, l) = z -  weighted_least_square(evalJ(z,l), G(z), wts)

    zk = copy(z0); zk1 = update(zk, l)
    deltaold = norm(zk1 - zk,2); zk = zk1
    
    cvg = false
    for ctr in 1:maxsteps
        zk1 = update(zk, l)
        delta = norm(zk1 - zk, 2)

        if delta > deltaold
            println("Growing delta. Best guess is being returned.")
            break
        end

        ## add extra abs(delta) < 100*eps() condition
        if delta^2 / (deltaold - delta) < tol || abs(delta) < 100*eps()
            cvg = true
            break
        end

        deltaold = delta
        zk=zk1
    end

    if !cvg println(""" 
Returning the initial estimates, as the
algorithm failed to improve estimates for the roots on the given
pejorative manifold.  
""") 
        return(z0) 
    end 
    return(zk1)
end


## Main interface to finding roots of polynomials with multiplicities
##
## The `multroot` function returns the roots and their multiplicities
## for `Poly` objects. It performs better than `roots` if the
## polynomial has multiplicities. 
##
## julia> x = poly([0.0]);
## julia> p = (x-1)^4 * (x-2)^3 * (x-3)^3 * (x-4)l
## julia> multroot(p)
## ([1.0,2.0,3.0,4.0],[4,3,3,1])
## ## For "prettier" printing, results can be coerced to a dict
## julia> [k => v for (k,v) in zip(multroot(p)...)]
## Dict{Any,Int64} with 4 entries:
##   1.0000000000000007 => 4
##   3.000000000000018  => 3
##   1.9999999999999913 => 3
##   3.999999999999969  => 1
## ## Large order polynomials prove difficult. We can't match the claims in Zeng's paper
## ## as we don't get the pejorative manifold structure right.
## julia> p = poly([1.0:10.0]);
## julia> multroot(p) ## should be 1,2,3,4,...,10 all with multplicity 1, but
## ([1.0068,2.14161,3.63283,5.42561,7.25056,8.81228,9.98925],[1,2,1,2,2,1,1])
##
## nearby roots can be an issue
## julia> delta = 0.0001  ## delta = 0.001 works as desired.
## julia> p = (x-1 - delta)*(x-1)*(x-1 + delta)
## julia> multroot(p)
## ([0.999885,1.00006],[1,2])
function multroot(p::Poly;
                  θ::Real=1.0,   # singular threshold, 1.0 is replaced by normf(A)*eps()
                  ρ::Real=1e-10, # initial residual tolerance
                  ϕ::Real=1e2,   # residual tolerance growth factor
                  δ::Real=1e-8   # passed to solve y sigma

                  )

    degree(p) == 0 && error("Degree of `p` must be atleast 1")
    
    if degree(p) == 1
        return(roots(p), [1])
    end

    ## if degree(p) == 2
    ##     a,b,c = coeffs(p)
    ##     discr = b^2 - 4a*c
    ##     if discr < 0
    ##         discr = Complex(discr, 0)
    ##     end
    ##     return  ( -2c / (-b - sqrt(discr)), -2c/(-b + sqrt(discr)))
    ## end
    
    p = Poly(float(coeffs(p)))  # floats, not Int

    u_j, v_j, w_j, residual= agcd(p, polyder(p),  ρ = ρ) 
    ρ = max(ρ, ϕ * residual)

    ## bookkeeping
    zs = roots(v_j)
    ls = ones(Int, length(zs))

    p0 = u_j

    while degree(p0) > 0
        if degree(p0) == 1
            z = roots(p0)[1]
            @compat _, ind = findmin(abs.(zs .- z))
            ls[ind] = ls[ind] + 1
            break
        end

        u_j, v_j, w_j, residual= agcd(p0, polyder(p0), ρ=ρ)

        ## need to worry about residual between
        ## u0 * v0 - monic(p0) and u0 * w0 - monic(Polynomials.polyder(p0))
        ## resiudal tolerance grows with m, here it depends on 
        ## initial value and previous residual times a growth tolerance, ϕ
        ρ = max(ρ, ϕ * residual)

        ## update multiplicities
        for z in roots(v_j)
            @compat _, ind = findmin(abs.(zs .- z))
            ls[ind] = ls[ind] + 1
        end

        ## rename
        p0 = u_j
    end


    if maximum(ls) == 1
        return(zs, ls)
    else
        zs = pejroot(p, zs, ls)
        return(zs, ls)
    end
end 

## Different interfaces

## can pass in vector too
multroot{T <: Real}(p::Vector{T}; kwargs...) = multroot(Poly(p); kwargs...)

## Can pass in function
function multroot(f::Function; kwargs...)
    p = Poly([0.0])
    try
        p = convert(Poly{Float64}, f)
    catch err
        error("The function does not compute a univariate polynomial")
    end
    multroot(p; kwargs...)        



    
end

## add funciton interface to Polynomials.roots
Polynomials.roots(f::Function) = roots(convert(Poly{Float64}, f))

