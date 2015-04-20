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

## poly p[] is a_0 + a_1 x^1 + a_2x^2 + ... + a_nx^n
## coeffs gives [a_0, ..., a_n], rcoeffs = [a_n, a_n-1, ..., a_1, a_0]
coeffs(p::Poly) = [p[i] for i in 0:degree(p)]
## rcoeffs needed as indexing in paper is an, an-1, ...
rcoeffs(p::Poly) = reverse(coeffs(p))
Base.norm(p::Poly, args...) = norm(coeffs(p), args...)
normf(A::Matrix) = sqrt(sum(A.^2))

leading_term(p::Poly) = p[degree(p)]
monic(p::Poly) = p/leading_term(p)

function Base.convert{T<:Integer}(::Type{Poly{T}}, p::Poly{Rational{T}})
    l = reduce(lcm, [x.den for x in p.a])
    q = l * p
    Poly(T[x for x in q.a], p.var)
end

## convert Poly <--> function
Base.convert(::Type{Function}, p::Poly) = x -> Polynomials.polyval(p,x)
## convert a function to a polynomial with error if conversion is not possible
QQR = Union(Int, BigInt, Rational{Int}, Rational{BigInt}, Float64)
function Base.convert{T<:QQR}(::Type{Poly{T}}, f::Function)
    try
        x = poly(zeros(T,1))
        out = f(x)
        if !isa(out, Poly)
            out = Poly([out])   # maybe a constant
        end
        out
    catch e
        rethrow(e)
    end
end
function Base.convert(::Type{Poly}, f::Function)
    ## try integers first, then float
    for T in [BigInt, Int, Float64]
        try
            fn = convert(Poly{T}, f)
            return(fn)
        catch e
        end
    end
    DomainError()
end


*{T, S}(A::Array{T,2}, p::Poly{S}) = Poly(A * rcoeffs(p))



## map monic(p) to a point in C^n
## p = 1x^n + a1x^n-1 + ... + an_1 x + an -> (a1,a2,...,an)
function p2a(p::Poly) 
    p = monic(p)
    rcoeffs(p)[2:end]
end


function weighted_least_square(A::Array, b::Vector, w::Vector)
    ## solve weighted least squares problem W*(Ax - b) = 0
    W = diagm(w)
    (W * A) \ (W * b)
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
                 wts::Union(Vector, Nothing)=nothing, # weight vector
                 tol = 1e-8,
                 maxsteps = 100
                      )
    
    a = p2a(p) #rcoeffs(monic(p))[2:end] # an_1, an_2, ..., a2, a1, a0
    
    if isa(wts, Nothing)
        wts = map(u -> min(1, 1/abs(u)), a)
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


## GCD code

## Various matrices that are needed:

## Cauchy matrix C_k(p), defn 2.1
## will be size (n+k) x k
function cauchy_matrix(p::Poly, k::Integer)
    n = degree(p) + 1
    out = zeros(eltype(p), n + k - 1, k)

    for i in 1:k
        out[(1:n) + (i-1), i] = rcoeffs(p)
    end
    out
end

## p21
## will be of size (n + k + 1) x (2k + 1)
function sylvester_discriminant_matrix(p::Poly, k::Integer)
    a = cauchy_matrix(Polynomials.polyder(p), k+1)
    b = cauchy_matrix(p, k)
    [a b]
end


## Lemma 4.1, (25)
## Jacobian F(u,v,w) = [p,p'] is J(u,v,w)
function JF(u::Poly, v::Poly, w::Poly)
    m, k = degree(u), degree(v)
    n = m + k
    a = cauchy_matrix(v,m+1)
    b = cauchy_matrix(w,m+1)

    c = cauchy_matrix(u,k+1)
    d = zeros(eltype(a), m+k, k+1)

    e = zeros(eltype(a), m+k+1, k)
    f = cauchy_matrix(u, k)

    [a c e; 
     b d f]
end

## \hat{S}_j(p) from page 21
function Shat(p::Poly, k::Integer)
    ## intersperse columns
    S = sylvester_discriminant_matrix(p, k)
    Sh = similar(S)
    for j in 1:(k+1)
        Sh[:, 1 + (j-1)*2] = S[:,j]
    end
    for j in 1:k
        Sh[:, j*2] = S[:, k + 1 + j]
    end
    Sh
end

## update Shat_k+1 from Shat_k, cf comment on p21
function Shatk1(Shatk::Array)
    m,n = size(Shatk)
    k = (n-1)/2                 # n = 2k + 1
    Sh = zeros(eltype(Shatk), m +1, n+2)
    Sh[1:m, 1:n] = Shatk
    Sh[:, n + 1] = vcat(0, Shatk[:,n-1])
    Sh[:, n + 2] = vcat(0, Shatk[:,n])
    Sh
end


## converge on right singular vector of A associated with singular value sigma (not a double value)
function lemma24(A::Matrix; tol::Real=1e-8)
    MAXSTEPS = 100
    ## how to @assert that sigma_2 > sigma_1?
    q,r = Base.qr(A)

    ## won't work if r is rank deficient
    if rank(r) < size(r)[2]
        λs, vs = eig(r)
        _, ind = findmin(abs(λs))
        return(λs[ind], vs[:,ind])
    end


    x = A'[:,1]                 # initial guess. Must not have Ax=0
    σ, σ1 = 1e8, 0

    function update(x)
        y = conj(r)' \ x
        z = r \ y
        x = z/norm(z,2)
        sigma = norm(r * x, 2)  # y/norm(z,2)
        (x, minimum(abs(sigma)))
    end

    ctr = 1


    ## how long do we update? Size is only issue, not accuracy so we iterate until we stop changing much
    ## ??? This needs engineering
    while (ctr < MAXSTEPS)
        x, σ1 = update(x)
        if abs((σ - σ1) / σ1) < 1.1
            break
        else
            σ = σ1
        end

        ctr += 1
    end

    if ctr < MAXSTEPS
        return(σ1, x)
    else
        error("no convergence?")
    end
end


## p21 implementation
## take p, return *guesses* for m, u, v, w 
## with: 
## m = number of distinct roots
## u * v = p
## u * w = monic(p')
## Basic idea is that m = n-k if S_k(p) is the first singular Sylvester Discriminant Matrix
## The computation solve_y_sigma(Sh) finds σ and y, with σ the smallest eigenvalue of "R" in Q,R.
## the initial guesses for u, v, w are derived from the eigenvalue y
##
## Main issue is the determination of what is 0. Paper suggests using θ * norm(p,2) but that
## fails to work for W_n (Wilcoxon poly) by n=5. 
function initial_gcd_system(p::Poly;
                    θ::Real=1e-12, # used to determine if singular value is essentialy 0
                    δ::Real=1e-12  # used to control solving of singular value
                    )

    j = 1
    n, normp = degree(p), norm(p, 2)
    Sh = Shat(p, j)
    σ, y = lemma24(Sh, tol = δ) 

    while (σ >  θ * normp) && (size(Sh)[2] <= size(Sh)[1])
        σ, y = lemma24(Sh, tol = δ)
        j += 1
        Sh = Shatk1(Sh)
    end


    if σ <  θ * normp

        ##  initial guess for u, v, w.from Section 4.1.3
        v = monic(Poly(reverse(y[1:2:length(y)])))
        w = monic(Poly(reverse(y[2:2:length(y)])))
            
        ## Solve C_(m+1)(v0) u0 = p for u0
        B = cauchy_matrix(v, degree(p) - degree(v) + 1)  
        wts = Float64[min(1, 1/abs(ak)) for ak in rcoeffs(p)]
    
        u = monic(Poly(reverse(weighted_least_square(B, rcoeffs(p), wts))))
    else
        ## no reduction
        u = Poly([1.0])
        v = p
        w = Polynomials.polyder(p)
    end

    return(degree(v), u, v, w)
end

## Algorithm 2

## compute approximate gcd starting with guess
## This is section 4.1.2
## we use agcd(p) below that calls this.
function agcd(p::Poly, q::Poly, u0::Poly, v0::Poly, w0::Poly;
              ρ::Real = 1e-10   # residual tolerance
)
    ## for monic polynomials p and q return
    ## p(x) = u(x) v(x)
    ## q(x) = u(x) w(x) with v, w co-prime

    ## should have monic(p) - un * vn close to 0
    ##             monic(q) - un * wn close to 0
    ## err residual on p25 (32), weighted L2 norm of [u0*v0-p, u0*w0-monic(p')]

    ## Return approximate values of u, v, w, err
    
    u, v, w = map(copy, (u0, v0, w0))
    ## not general, we use q here
    q = Roots.monic(Polynomials.polyder(p))
    m,n,k = map(Polynomials.degree, (u, v, w))
    wts =  map(pj -> 1/max(1, abs(pj)), vcat(rcoeffs(p), rcoeffs(q)))
       
    ## compute F(u,v,w) - [p, p'] = [u*v, u*w] - [p, p']
    Fmp(p,q,u,v,w) = [rcoeffs(u*v); rcoeffs(u*w)] - [rcoeffs(p); rcoeffs(q)]
    residual_error(p,q,u,v,w) = norm(Fmp(p,q,u,v,w) .* wts, Inf)
        
    ## iterate to solve
    ## uvw_j+1 = uvw_j - Jw[F - p]

    function update(u, v, w)
        A = Roots.JF(u, v, w)
        b = Fmp(p,q,u,v,w)

        inc = Roots.weighted_least_square(A, b, wts)

        x = vcat(rcoeffs(u), rcoeffs(v), rcoeffs(w))
        x = x - inc

        m,n = Roots.degree(u), Roots.degree(v)
        u = Poly(reverse(x[1:(1+m)]))
        v = Poly(reverse(x[(m+2):(m+n+2)]))
        w = Poly(reverse(x[(m+2+n+1):end]))
        
        err = residual_error(p,q,u,v,w)
        (err, u, v, w)
    end

    err0, err1 = Inf, residual_error(p,q,u,v,w)

    ctr = 1
    while true
        try
            err1, u1, v1, w1 = update(u, v, w)
            if err1 < err0
                err0, u, v, w = err1, u1, v1, w1
            else
                break
            end
        catch e
            break               # possibly singluar
        end

        (ctr += 1) > 20 && break
    end

    if err0 < ρ
        u,v,w = map(Roots.monic, (u,v,w))
        (u,v,w,err0)
    else
        ## failed to converge, so we return the initial guess
        (u0,v0,w0, err0)
    end

end

    


## Preconditioning code
function geometric_mean(as)
    as = abs(as)
    as = filter(x -> x >= 1e-14, as)
    prod(as)^(1/length(as))
end

function ratio(p,q)
    is_nonzero(x) = x >= 1e-14
    as = abs(filter(is_nonzero, p.a))
    bs = abs(filter(is_nonzero, q.a))

    max(maximum(as), maximum(bs)) / min(minimum(as), minimum(bs))
end

function precondition(p,q)
    ## find alpha, gamma that minimize ratio of max coefficients to min coefficients, for getting zeros
    alphas = [2.0^i for i in -5:5]
    phis = [2.0^i for i in -5:5]
    out = (1,1)
    m = ratio(p,q)

    for α in alphas
        for ϕ in phis
            r = ratio(polyval(p, ϕ * variable(p)), α * polyval(q, ϕ * variable(q)))
            if r < m
                out = (α, ϕ)
            end
        end
    end

    α, ϕ = out
    p = Roots.polyval(p, ϕ * Roots.variable(p))
    q = α * Roots.polyval(q, ϕ * Roots.variable(q))

    p = p * (1 / geometric_mean(p.a))
    q = q * (1 / geometric_mean(q.a))

    p,q, ϕ, α
    
end

## return (u,v,w,residual) for gcd of (p, p') (u*v=p, u*w=p')
## We precondition here ala https://who.rocq.inria.fr/Jan.Elias/pdf/je_masterthesis.pdf
## (We should write a general agcd(p,q) function, as this is specific to the task of agcd(p,p'))
function agcd(p::Poly;
              δ::Real=1e-8,
              ρ::Real=1e-10)
    q = Polynomials.polyder(p)
    pp, qq, Φ, α = precondition(p, q)
    θ = normf(Shat(pp,1)) * eps()
    m, u_j0, v_j0, w_j0 = initial_gcd_system(pp, θ = θ, δ = δ)
    u_j, v_j, w_j, residual= agcd(pp, qq, u_j0, v_j0, w_j0, ρ = ρ) 

    x = (1/Φ) * variable(p)
    (polyval(u_j, x), polyval(v_j, x), polyval(w_j, x), residual)
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

    
    p = Poly(float(coeffs(p)))  # floats, not Int
    q = Polynomials.polyder(p)

    ## initial
    ##    m, u_j0, v_j0, w_j0 = initial_gcd_system(p, θ = θ, δ = δ)
    ##    u_j, v_j, w_j, residual= agcd(p, q, u_j0, v_j0, w_j0, ρ = ρ)
    u_j, v_j, w_j, residual= agcd(p, ρ = ρ) 
    ρ = max(ρ, ϕ * residual)

    ## bookkeeping
    zs = roots(v_j)
    ls = ones(Int, length(zs))

    p0 = u_j

    while degree(p0) > 0
        if degree(p0) == 1
            z = roots(p0)[1]
            _, ind = findmin(abs(zs .- z))
            ls[ind] = ls[ind] + 1
            break
        end

        ##m, u_j0, v_j0, w_j0 = initial_gcd_system(p0, θ = θ, δ = δ)
        ##u_j, v_j, w_j, residual= agcd(p0, Polynomials.polyder(p0), u_j0, v_j0, w_j0, ρ = ρ) 
        u_j, v_j, w_j, residual= agcd(p0, ρ=ρ)
        ## need to worry about residual between
        ## u0 * v0 - monic(p0) and u0 * w0 - monic(Polynomials.polyder(p0))
        ## resiudal tolerance grows with m, here it depends on 
        ## initial value and previous residual times a growth tolerance, ϕ
        ρ = max(ρ, ϕ * residual)

        ## update multiplicities
        for z in roots(v_j)
            _, ind = findmin(abs(zs .- z))
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
    try
        p = convert(Poly{Float64}, f)
        multroot(p; kwargs...)
    catch e
        error("The function does not compute a univariate polynomial")
    end
end

## add funciton interface to Polynomials.roots
Polynomials.roots(f::Function) = roots(convert(Poly{Float64}, f))

