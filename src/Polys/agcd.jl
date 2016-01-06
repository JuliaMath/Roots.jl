##
## Find an approximate gcd (agcd) for two polynomials
## The `gcd` function can exactly handle Poly{Int} if Rational{BigInt} is used:
## p = poly([1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,3,3])
## q = Poly(map(u->big(u)//1, coeffs(p)), p.var)
## gcd(p, polyder(p)) # Poly(5.960464477539063e-8)
## gcd(q, polyder(q)) # degree 17, as expected
## But what to do with Poly{Float}? There floating point issues will be a problem,
## This create `agcd` for approximate gcd:
## u,v,w,err = agcd(p, polyder(p))
## degree(u)  # 17
## roots(v)   # 1,2,3


## reimplement Zeng's gcd code for approximate gcds

## Special matrices
function cauchy_matrix{T}(p::Poly{T}, k::Integer)
    n = degree(p) + 1
    out = zeros(T, n + k - 1, k)
    for i in 1:k
        out[(1:n) + (i-1), i] = rcoeffs(p)
    end
    out
end


function sylvester_matrix(p::Poly, q::Poly; k::Int=0)
    @assert k >= 0
    n,m = degree(p), degree(q)
    if n < m
        p,q = q,p
        n,m = m,n
    end

    del = n - m
    i,j = k + del, k
    if k == 0
        i,j = n,m
    end
    hcat(cauchy_matrix(q, i), cauchy_matrix(p, j))
end


## preconditioning code
## taken from https://who.rocq.inria.fr/Jan.Elias/pdf/je_masterthesis.pdf
function geometric_mean(a::Vector, epsilon=Base.eps())
    a = filter(x -> abs(x) > epsilon, a)
    n = length(a)
    prod(abs(a) .^ (1/n))
end

function ratio(p,q, atol=Base.eps(), rtol=Base.eps())
    is_nonzero(x) = !isapprox(x, 0; rtol=rtol, atol=atol)
    as = abs(filter(is_nonzero, p.a))
    length(as) == 0 && return Inf

    bs = abs(filter(is_nonzero, q.a))
    length(bs) == 0 && return Inf

    max(maximum(as), maximum(bs)) / min(minimum(as), minimum(bs))
end

## find alpha, gamma that minimize ratio of max coefficients to min coefficients, for getting zeros
## 1.12 writes this as a linear programming problem, we just ...
function precondition(p::Poly,q::Poly)
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
    
    p = polyval(p, ϕ * variable(p))
    q = α * polyval(q, ϕ * variable(q))

    p = p * (1/geometric_mean(coeffs(p)))
    q = q * (1/geometric_mean(coeffs(q)))
    
    p,q, ϕ, α
    
end


## converge on right singular vector of A associated with singular value sigma (not a double value)
## we return if sigma < tol; delta \appro 0;
## how to @assert that sigma_2 > sigma_1?
function lemma24(A::Matrix; θ::Real=1e-8)
    const MAXSTEPS = 100

    Q,R = Base.qr(A)
    if rank(R) < size(R)[2]
        λs, vs = eig(R)
        _, ind = findmin(abs(λs))
        return(λs[ind], vs[:,ind])
    end

    x = rand(size(A)[2]); x/norm(x,2)   # use random initial guess
    σ, σ1 = 1e8, Inf

    function update(x)
        y = conj(R)' \ x
        z = R \ y
        x = z/norm(z,2)
        sigma = norm(R * x, 2)  # y/norm(z,2)
        (x, minimum(abs(sigma)))
    end

    ## how long do we update? Size is only issue, not accuracy so we iterate until we stop changing much
    for _ in 1:MAXSTEPS
        x, σ1 = update(x)
        ((σ1 < θ) | (abs((σ - σ1) / σ1) < 1.1)) && break
        σ = σ1
    end

    return(σ1, x)
end


"""
Return k, u,v,w where k reveals rank; u*v \approx p; u*w \approx q; v & w coprime

Following Zeng, section 4, this could be made more efficient.
In Zeng, theta=1e-8. Here agcd uses 1e-12, which is an improvement on some...
"""
function reveal_rank(p::Poly, q::Poly, theta=1e-8)
    n,m = map(degree, (p, q))
    if n < m
        p,q = q,p
        n,m=m,n
    end
    
    for k in 1:m
        A = sylvester_matrix(p,q,k=k)
        psi, x = lemma24(A)
        if abs(psi) < norm(p,2) * theta #  norm(A,2)?
            ## degree u = n - k; degree(v) = k
            v = monic(Poly(reverse(x[1:(length(x)-k)]), p.var))
            w = monic(Poly(reverse(x[(length(x)-k+1):end]), p.var))
            u = monic(Poly(reverse(cauchy_matrix(v, degree(p) - degree(v) + 1) \ rcoeffs(p)), p.var))
            return (k, u, v,w)
        end
    end
    return (n+m, Poly(ones(eltype(coeffs(p)),1)), p, q)
end    


## Lemma 4.1, (25)
## solve weighted least squares problem W*(Ax - b) = 0
function weighted_least_square(A, b, w)
    W = diagm(w)
    (W * A) \ (W * b)
end

## Jacobian F(u,v,w) = [p,p'] is J(u,v,w)
function JF{T}(u::Poly{T}, v::Poly, w::Poly)
    j, k, l= degree(u), degree(v), degree(w)
    n, m = j + k, j + l

    a = cauchy_matrix(v, n+1-k)
    b = cauchy_matrix(u, n+1-j)
    c = zeros(T, n+1, m+1-j)

    d = cauchy_matrix(w, m+1-l)
    e = zeros(T, m+1, n+1-j)
    f = cauchy_matrix(u, m+1-j)

    
    A = hcat(a,b,c)
    B = hcat(d,e,f)
    vcat(A,B)
end
## compute F(u,v,w) - [p, p'] = [u*v, u*w] - [p, p']
Fmp(p,q,u,v,w) = [rcoeffs(u*v); rcoeffs(u*w)] - [rcoeffs(p); rcoeffs(q)]
## error in estimate for p=u*v,q=u*w for some weights
residual_error(p,q,u,v,w, wts=ones(degree(p) + degree(q) + 2)) = norm(Fmp(p,q,u,v,w) .* wts, 2)
function agcd_update(p, q, u, v, w, wts)
    m,n = map(degree, (u,v))
    A = JF(u, v, w)
    b = Fmp(p,q,u,v,w)
    inc = weighted_least_square(A, b, wts)
    
    x = vcat(rcoeffs(u), rcoeffs(v), rcoeffs(w))
    x = x - inc
    
    u = Poly(reverse(x[1:(1+m)]), p.var)
    v = Poly(reverse(x[(m+2):(m+n+2)]), p.var)
    w = Poly(reverse(x[(m+2+n+1):end]), p.var)
    
    err = residual_error(p,q,u,v,w, wts)
    (u, v, w, err)
end

"""


Find an approximate GCD for polynomials `p` and `q` using an algorithm of [Zeng](http://www.ams.org/journals/mcom/2005-74-250/S0025-5718-04-01692-8/home.html).


Returns u,v,w, err where:

* `u*v \approx monic(p)`
* `u*w \approx monic(q)`
* The total residual error in these approximations is bounded by `err`.

Further, 
* `v` and `w` should share no common roots (`u` is a gcd of `u*v` and `u*w`)
* `roots(v)` should exhaust unique values of `roots(p)`.

The tolerances are:

* theta: passed to the reveal_rank function. In Zeng 1e-8 is used. Here 1e-12 seems to work better?
* \rho: if the residual error does not come below this mark, then we use the initial guess


"""
function agcd{T,S}(p::Poly{T}, q::Poly{S}=p';
              theta = 1e-12,    # reveal_rank tolerance. (1e-8 in paper, this seems better?)
              ρ::Real = 1e-10   # residual tolerance
              )

    n, m = map(degree,(p,q))
    if m > n
        p,q=q,p
    end

    if m == 0
        return (Poly(ones(1)), p, q, 0)
    end
    
    p0,q0 = map(copy, (p,q))
    p,q, phi, alpha = precondition(p,q) 

    k, u, v, w = reveal_rank(p, q, theta)
    u0,v0,w0 = map(copy, (u,v,w))
    
    m,n,k = map(degree, (u, v, w))
    wts =  map(pj -> 1/max(1, abs(pj)), vcat(rcoeffs(p), rcoeffs(q)))
       
    ## iterate to solve
    ## uvw_j+1 = uvw_j - Jw[F - p]
  
    err0, err1 = Inf, residual_error(p,q,u,v,w, wts)
    for ctr in 1:20
        try
            u1, v1, w1, err1 = agcd_update(p, q, u, v, w, wts)
            if err1 < err0
                err0, u, v, w = err1, u1, v1, w1
            else
                break
            end
        catch err
            break                       # possibly singular
        end
    end

    if err0 > norm(sylvester_matrix(p,q), 2) * ρ
        u,v,w = map(monic, (u0, v0, w0))  ## failed to converge, so we return the initial guess
    end

    x = (1/phi) * variable(p) # reverse preconditioning
    u,v,w = map(monic, (polyval(u, x), polyval(v, x), polyval(w, x)))

    (u,v,w, residual_error(p0,q0,u,v,w))
end

    

    

