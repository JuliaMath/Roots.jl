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

using Polynomial
import Base.*

## poly funs to add to Polynomial.jl
polyder(p::Poly) = polydir(p)
degree(p::Poly) = length(p) - 1
coeffs(p::Poly) = p.a
Base.norm(p::Poly, args...) = norm(p.a, args...)
leading_term(p::Poly) = p.a[1]
monic(p::Poly) = p/leading_term(p)
monic!(p::Poly) = (p.a = p.a/leading_term(p))
Base.start(p::Poly) = 1
Base.next(p::Poly, i) = (p[i], i + 1)
Base.done(p::Poly, i) = i > length(p)
Base.convert(::Type{Poly{Float64}},p=Poly{Int64}) = Poly(float(p.a))
Base.convert(::Type{Function}, p::Poly) = x -> Polynomial.polyval(p,x)
*{T, S}(A::Array{T,2}, p::Poly{S}) = Poly(A * p.a)




function weighted_least_square(A::Array, b::Vector, w::Vector)
    ## solve weighted least squares problem W*(Ax - b) = 0
    W = diagm(w)
    (W * A) \ (W * b)
end

## get value of gl(z). From p16
function evalG(z::Vector, l::Vector)
    m = length(z)

    length(z) == length(l) || throw("Length mismatch")

    s = Poly([1.0])
    for i in 1:m
        s = s * Poly([1.0,-z[i]])^l[i]
    end
    s.a[2:end]
end

## get jacobian J_l(z), p16
function evalJ(z::Vector, l::Vector)
    m = length(z) 

    m == length(l) || throw("Length mismatch")

    u = prod([Poly([1.0, -z[j]])^(l[j]-1) for j in 1:m]) ## Pi (1-z)^(l-1)

    J = zeros(eltype(z), sum(l), m)
    for j in 1:m
        s = -l[j] * u
        for i in 1:m
            if i != j
                s = s * Poly([1.0, -z[i]])
            end
        end
        J[:,j] = s.a
    end
    J
end

## Gauss Newton iteration to solve weighted least squares problem
## G_l(z) = a, where a is related to monic version of polynomial p
## l is known multiplicity structure of polynomial p = (x-z1)^l1 * (x-z2)^l2 * ... * (x-zn)^ln
## Algorithm I, p17
function pejroot{T<:Int}(p::Poly, z0::Vector, l::Vector{T};
                 w::Union(Vector, Nothing)=nothing, # weight vector
                 tol = 1e-8,
                 maxsteps = 100
                      )
    
    a = p.a[2:end]/p.a[1]       # monic
    
    if isa(w, Nothing)
        w = map(u -> min(1, 1/abs(u)), a)
    end
    W = diagm(w)

    ## Could make this more efficient I'm sure
    G(z) = (evalG(z, l) - a)
    function Jplus(z) 
        Jk = evalJ(z, l)
        B = conj(Jk)' * W^2 
        (B * Jk) \ B
    end
    update(z, l) = z - Jplus(z) * G(z)

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

        if delta^2 / (deltaold - delta) < tol || abs(delta) < 100*eps()
            cvg = true
            break
        end

        deltaold = delta
        zk=zk1
    end

    if !cvg
        println("Failed to converge in $maxsteps steps. Returning last guess.")
    end
    zk1
end


## GCD code

## Various matrices that are needed:

## Cauchy matric C_k(p), defn 2.1
function cauchy_matrix(p::Poly, k::Integer)
    n = degree(p) + 1
    out = zeros(eltype(p.a), n + k - 1, k)

    for i in 1:k
        out[(1:n) + (i-1), i] = p.a
    end
    out
end

## p21
function sylvester_discriminant_matrix(p::Poly, k::Integer)
    a = cauchy_matrix(polydir(p), k+1)
    b = cauchy_matrix(p, k)
    [a b]
end


## lemma 4.1, (25)
## Jacobian of F(u,v,w) -> [p, p;]
function JF(u::Poly, v::Poly, w::Poly)
    m, k = degree(u), degree(v)
    n = m + k
    
    a = cauchy_matrix(v,m)
    b = cauchy_matrix(w,m)
    c = cauchy_matrix(u,k)
    d = zeros(eltype(a), n-1, k)
    e = zeros(eltype(a), n, k-1)
    f = cauchy_matrix(u, k-1)

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

## Algorithm 2

## compute approximate gcd starting with guess
function agcd(p::Poly, q::Poly, u0::Poly, v0::Poly, w0::Poly)
    ## for monic polynomials p and q return
    ## p(x) = u(x) v(x)
    ## q(x) = u(x) w(x) with v, w co-prime

    ## should have monic(p) - un * vn close to 0
    ##             monic(q) - un * wn close to 0
    ## err residual on p25 (32), weighted L2 norm of [u0*v0-p, u0*w0-monic(p')]

    ## Return approximate values of u, v, w, err
    
    ## get degrees and make monic, compute weights
    m,n,k = map(degree, (u0, v0, w0))
    p,q,u0,v0,w0 = map(monic,  [p,q, u0, v0, w0])
    wts = map(pj -> 1/max(1, abs(pj)), vcat(p.a[2:end], q.a[2:end]))

    s = u0 * v0 - p
    t = u0 * w0 - q

    b = vcat(s.a[2:end], t.a[2:end])
    xk = vcat(u0.a[2:end], v0.a[2:end], w0.a[2:end])
    max_coeff = [norm(b.*wts, Inf)]
    
    j = 1
    ## refine values
    while j > 0
        A = JF(u0, v0, w0)
        xk1 = xk - weighted_least_square(A, b, wts)

        u = Poly([1.0, xk1[1:m]]) 
        v = Poly([1.0, xk1[(m+1):(m+n)]]) 
        w = Poly([1.0, xk1[(m+n+1):(m+n+k)]]) 

        s = u*v - p
        t = u*w - q

        b = vcat(s.a[2:end], t.a[2:end])
        push!(max_coeff, norm(b.*wts, Inf))
        j += 1

        if max_coeff[j] >= max_coeff[j-1]
            return(u, v, w, pop!(max_coeff))
        end

        u0 = u; v0=v; w0=w
        xk = xk1
    end
end

## Helper function to gcd_degree
## return sigma=smallest eigen value and eigen vector of A^H*A
## lemma 2.4, p5
function solve_y_sigma(A::Array; tol=1e-8, maxsteps::Int=200)
    q, r= Base.qr(A)
    if abs(det(r)) < 1e2 * eps()      # basically 0
        ## Lemma 2.4 does not apply. Return 0 and eigenvalue for A^H*A
        B = conj(A)' * A
        (vals, vecs) = eig(B)
        (val, i) = findmin(abs(vals))
        return (val, vecs[:,i])
    end


    xk = rand(size(A)[2]) * 2 - 1 # in (-1,1)

    function update(x)
        y = conj(r)' \ x
        z = r \ y
        x = z / norm(z,2)
        sigma = norm(r*x,2)
        (x, sigma)
    end

    (xk, sigmak) = update(xk)
    (xk, sigmak1) = update(xk)
    cvg = false

    for ctr = 1:maxsteps
        xk, sigmak1 = update(xk)
        if abs(sigmak1 - sigmak) < tol
            cvg = true
            break
        end
        sigmak = sigmak1
    end
    if !cvg
        println("sigmak=", sigmak, "y=", xk)
        throw("Could not solve for y and sigma")
    end
    (sigmak1, xk)
end

## compute degree of gcd and return improved guesses u0, v0, w0 for p, p'
## used in 26 to get d=degree(v0), roots(v0)
## p21
## use implicit inverse iteration (2.4) to find smallest singular value
## zetaj of Shat(p,j) and corresponding right singluar vector yj
function gcd_degree(p::Poly; 
                    theta::Real=1e-8, # for checking if sigma small
                    delta::Real=1e-8, # passed in iteration to find sigma solve_y_sigma
                    rho::Real=1e-10, phi::Real=1e2 # to refine um, vm,wm
                    )
    if degree(p) <= 1
        return (degree(p), Poly([1.]), p, monic(polyder(p)))
    end

    normp = norm(p)
    Sh = Shat(p, 1)
    
    residual, sigma = Inf, Inf
    l = 1

    while l > 0
        sigma, y = solve_y_sigma(Sh, tol=delta) 
        if sigma < theta * normp
            ## compute um,vm,wm ...
            ## odd entries form v0, even entries form w0
            v0 = monic(Poly(y[filter(isodd, [1:length(y)])]))
            w0 = monic(Poly(y[filter(iseven, [1:length(y)])]))
            
            ## Solve C_(m+1)(v0) u0 = p for u0
            B = cauchy_matrix(v0, degree(p) - degree(v0) + 1)  
            wts = [min(1, 1/abs(ak)) for ak in p]
            
            u0 = monic(Poly(weighted_least_square(B, p.a, wts)))
            ## u0, v0, w0 are initial guesses, these are now improved
            ## with agcd
            um, vm, wm, residual= agcd(p, polyder(p), u0, v0, w0) # refine tolerance?

            if residual < rho * normp
                ## return degree, gcd factorization u,v,w
                return (degree(um), um, vm, wm)
            else
                rho = phi*rho
            end
        end

        Sh = Shatk1(Sh)
        l = l + 1

        if size(Sh)[2] > size(Sh)[1]
            ## won't have a solution -- more columns than rows (DimensionMismatch("Matrix cannot have less rows than columns"))
            ## return m=1, u=1, v=p and w = p'
            return (0, Poly([1.0]), p, polydir(p))
        end
        
    end
end

## helper function to disambiguate roots and multiplicities
## find roots and multiplicites for roots
function find_fuzzy(zs)
    rs = zs[1]
    os = vcat(zs[2:end]...)
    l = ones(Int, length(rs))

    for r in os
        a, ind = findmin(abs(rs - r))
        l[ind] += 1
    end
    return (rs, l)
end

## Main interface to finding roots of polynomials with multiplicities
##
## The `multroot` function returns the roots and their multiplicities
## for `Poly` objects. It performs better than `roots` if the
## polynomial has multiplicities. 
##
## Call as
##   julia> p = Poly([1, -7, 17, -17, 6])  # (x-1)^3 * (x-2) * (x-3)
##   julia> z, l = multroot(p) ## z=[0.9999999999999994,2.0000000000000058,2.9999999999999947], l=[2,1,1]
##   julia> roots(p)           ## [0.9999999491767974, 1.0000000508232112, 1.9999999999999756, 3.0000000000000187]
##   ## About the same. Now for something more challenging  ## (x-1)^2*(x-2)^2*(x-3)^4
##   julia> p = Poly([1.0, -18.0, 139.0, -600.0, 1579.0, -2586.0, 2565.0, -1404.0, 324.0])
##   julia> roots(p)             
##            1.0-0.0im
##            1.0-0.0im
##            2.0-0.0im
##            2.0-0.0im
##   2.99828-0.00172297im
##   2.99828+0.00172297im
##   3.00172-0.00172259im
##   3.00172+0.00172259im
##   julia> z,l = multroot(p) ## ([1.0000000000000004,1.999999999999997,3.0000000000000013],[2,2,4])
##
##  Answers have slight randomness possible, as one algorithm has random initial starting guess
##  TODO: Not sure I have tolerances working properly. 
##  Somethings don't work well:
##  julia> p = Poly([1,2,1]) # (x-1)^2
##  julia> multroot(p^14)    # ([-1.0],[28])
##  julia> multroot(p^15)    # ([1.1273742514712568,-1.1326950594766867],[18,12]) # failed to converge
##  ## though still better than `roots` in some sense: 
##  julia> using Winston
##  julia> r = roots(p^15); plot(real(r), imag(r)) ## way off
##
function multroot(p::Poly;
                 theta::Real=1e-8,  # singular threshold
                 rho::Real = 1e-10, # residual tolerance
                 phi::Real = 1e2          # residual growth factor
                 )

    p = Poly(float(p.a)) 
    u = [p]
    d = Int[]
    zs = Any[]
    p0 = u[end]

    while true
        if degree(p0) == 0 
            break
        end
        if degree(p0) == 1
            push!(d, 1)
            push!(zs, roots(p0))
            break
        end
        (m, u0, v0, w0) = gcd_degree(p0; theta=theta)

        ## need to worry about residual between
        ## u0 * v0 - monic(p0) and u0 * w0 - monic(polyder(p0))
        push!(d, degree(v0))
        push!(zs, roots(v0))
        push!(u, u0)
        p0 = u[end]
    end

    ## now finish off. We have candidates for the roots zs, and multiplicities l
    ## if there are multiple roots we refine candidates through `pejroot`
    if length(d) == 1
        return(zs[1], ones(Int, length(zs[1])))
    else
        zs, l = find_fuzzy(zs)
        zs = pejroot(p, zs, l)
        return(zs, l)
    end
end        
  
## can pass in vector too
multroot{T <: Real}(p::Vector{T}; kwargs...) = multroot(Poly(p); kwargs...)