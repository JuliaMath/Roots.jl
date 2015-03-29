## Find real zeros of a polynomial
##
## We do this two ways, depending if an exact gcd can be found (Z[x] or Q[x])
##
## For polynomials over Z[x], Q[x], we can use a modification of
## [Vincent's theorem](http://en.wikipedia.org/wiki/Vincent's_theorem)
## taken from Efficient isolation of polynomial’s real roots by
## Fabrice Rouillier; Paul Zimmermann
##

###   Things for Polynomials   #####
## return "x"
variable(p::Poly) = poly(zeros(eltype(p),1), p.var)

## make bounds work for floating point or rational.
_iszero{T <: Union(Integer, Rational)}(b::T; kwargs...) = b == 0
_iszero{T<:FloatingPoint}(b::T; xtol=1) = abs(b) <= 2*xtol*eps(T)


## extend to evaluate with other polynomials
## Allows for composition of polys
import Polynomials: polyval
function polyval{T<:Number, S<:Number}(p::Poly{T}, u::Poly{S})
    ## mangle T, S???
    p.var == u.var || DomainError() # "Symbols must match"
    degree(p) == 0 && return(p)
    y = p[degree(p)]
    for i = (degree(p)-1):-1:0
        y = p[i] + u*y
    end
    return y
end

## return (q,r) with p(x) = (x-c)*q(x) + r using synthetic division
function Base.divrem(p::Poly, c::Number)
    ps = copy(p.a)                    # [p0, p1, ..., pn]
    qs = eltype(p)[pop!(ps)]           # take from right
    while length(ps) > 0
        unshift!(qs, c*qs[1] + pop!(ps))
    end
    r = shift!(qs)
    Poly(qs, p.var), r
end

## p(x) = q(x)*(x-c)^k for some q,k (possibly 0). Return maximal k and corresponding q
function multiplicity(p::Poly, c::Number)
    k = 0
    q, r = divrem(p,c)
    while _iszero(r)
        p = q
        q,r = divrem(p, c)
        k = k + 1
    end
    p, k
end

## Our Poly types for which we can find gcd
QQ = Union(Int, BigInt, Rational{Int}, Rational{BigInt})
BB = Union(BigInt, Rational{BigInt})

## Here computations are exact, as long as we return poly in Q[x]
function Base.divrem{T<:QQ, S<:QQ}(a::Poly{T}, b::Poly{S})
    degree(b) == 0 && (b[0] == 0 ? error("b is 0") : return (a/b[0], 0*b))
    
    lc(p::Poly) = p[degree(p)]

    a.var == b.var || DomainError() # "symbols must match"
    R = promote_type(T, S)
    x = poly(zeros(R,1), a.var)
    q, r= 0*a, a
    d,c = degree(b), lc(b)

    while degree(r) >= d
        s = lc(r)//c * x^(degree(r)- d)
        q, r = q + s, r - s*b
    end
    q,r
end

## gcd over Rational{BigInt}
function bgcd{T <: BB, S<: BB}(a::Poly{T}, b::Poly{S})
    (degree(b) == 0 && b[0] == 0) ? a : bgcd(b, divrem(a,b)[2])
end
##################################################
## Real zeros of p in Z[x] or Q[x], using VCA method
## cf. E$cient isolation of polynomial’s real roots by Fabrice Rouillier, Paul Zimmermann


## Polynomial transformations
" `R(p)` finds  `x^n p(1/x)` which is a reversal of coefficients "
R(p) = Poly(reverse(p.a), p.var)

" `p(λ x)`: scale x axis by λ "
Hλ(p, λ=1//2) = polyval(p, λ * variable(p))

" `p(x + λ)`: Translate polynomial left by λ "
Tλ(p, λ=1)   = polyval(p, variable(p) + λ)

" shift and scale p so that  [c,c+1]/2^k -> [0,1] "
function Pkc(p::Poly, k, c)
    n = degree(p)
    2^(k*n) * Hλ(Tλ(p, c/2^k), 1/2^k)
end


## Upper bound on size of real roots that is tighter than cauchy
## titan.princeton.edu/papers/claire/hertz-etal-99.ps
function upperbound(p::Poly)
    q, d = monic(p), degree(p)
    
    d == 0 && error("degree 0 is a constant")
    d == 1 && abs(q[0])


    a1 = abs(q[d-1])
    B = maximum([abs(q[i]) for i in 0:(d-2)])

    a,b,c = 1, -(1+a1), a1-B
    (-b + sqrt(b^2 - 4a*c))/2
end



"""

Use Descarte's rule of signs to compute number of *postive* roots of p

"""
¬(f::Function) = x -> !f(x)
function descartes_bound(p::Poly)
    bs = filter(¬_iszero, p.a)
    
    ## count terms
    ctr = 0
    bs = map(sign,bs)
    last = bs[1]
    for i in 2:length(bs)
        bs[i] == 0 && continue
        if last != bs[i]
            last = bs[i]
            ctr += 1
        end
    end
    ctr
end

## Descarte's upperbound on possible zeros in (0,1)
function DesBound(p::Poly)
    q = Tλ(R(p),1)
    descartes_bound(q)
end


## Interval with [c/2^k, (c+1)/2^k)
immutable Intervalkc
    k::Int
    c::Int
end

## convenience for two useful functions
DesBound(p,node::Intervalkc) = DesBound(Pkc(p,node)) 
function Pkc(p::Poly, node::Intervalkc)
    k,c = node.k, node.c
    Pkc(p, k, c)
end

## how we store the state of our algorithm to find zero in [0,1]
type State
    Internal                    # DesBound > 1
    Exact                       # c/2^k is a root
    Isol                        # DesBound == 1
    p::Poly
end
State(p) = State(Intervalkc[], Set{Rational}(), Intervalkc[], p)


## from naming convention of paper
function initTree(st::State)
    node = Intervalkc(0,0)
    ## check 0 and 1
    for x in [0//1, 1//1]
        if _iszero(polyval(st.p,x))
            push!(st.Exact, x)
            st.p,k = multiplicity(st.p, x)
        end
    end

    n = DesBound(Pkc(st.p, node))
    if n == 0
        ## nothing all done
    elseif n == 1
        push!(st.Isol, node)
    else
        push!(st.Internal, node)
    end
end

## get next node. Need to check that length(st.Internal) > 0, as it would throw error if no more Internal nodes
getNode(st::State) = shift!(st.Internal)

## add successors to [c,c+1]/2^k -> [2c, 2c+1]/2^(k+1) & [2(c+1), 2(c+1)+1]/2^(k+1) and check if c+1
function addSucc(st::State, node)
    k,c = node.k, node.c
    l = Intervalkc(k+1, 2c  )
    r = Intervalkc(k+1, 2c+1)
    p,m = multiplicity(st.p, (2c+1)//2^(k+1)); m > 0 && push!(st.Exact,  (2c+1)//2^(k+1))
    p,m = multiplicity(st.p, (2c+2)//2^(k+1)); m > 0 && push!(st.Exact,  (2c+2)//2^(k+1))
    for i in [l,r]
        n = DesBound(Pkc(p,i))
        n == 1 && push!(st.Isol, i)
        n > 1 && push!(st.Internal, i)
    end
end

## update state to resolve for roots in [0,1]
## st should have st.Interval = [], st.Exact and st.Isol
## This is basic search, not the more complicated one of the paper. To modify search
## strategy, the push!(st.Internal, i) would be changed in `addSucc`.
function find_in_01(p::Poly)
    st = State(p)
    
    initTree(st)
    for x in [0,1]
        p,k = multiplicity(p,x)
        k > 0 && push!(st.Exact, x)
    end
    while length(st.Internal) > 0
        node = getNode(st)
        node.k > 30 && error("too small?") # XXX arbitrary. Theory says a delta exists, this is δ > 1/2^30
        addSucc(st, node)
    end
    st
end

## get root from an Isol interval
function getIsolatedRoot(p::Poly, a::Real, b::Real)
    pa, pb = [polyval(p,x) for x in (a,b)]
    pab = sign(pa * pb)
    if pa == 0 p,k = multiplicity(p,a) end
    if pb == 0 p,k = multiplicity(p,b) end
    pa, pb = [polyval(p,x) for x in (a,b)]
    pa * pb >= 0 && println("DEBUG: $(sign(pa*pb)), $pab: Find root of $p in [$a, $b] will fail!!")
    pa * pb < 0 ? Roots.fzero(p, [a,b]) : -999
end

function getIsolatedRoot(p::Poly, i::Intervalkc)
    a,b = [i.c, i.c+1]//2^(i.k)
    c = (a+b)//2
    polyval(p,c) == 0 && return(c) # keeps rational if c/2^k type
    getIsolatedRoot(p, a, b)
end

## return positive roots of a square-free polynomial
function pos_roots(p::Poly)

    rts = Any[]

    ## in case 0 is a root, we remove from p. Should be unnecessary, but ...
    p,k = multiplicity(p, 0)

    M = Roots.upperbound(p)

    for k in 0:ifloor(M)  ## we consider [k, k+1) for roots until [k, ∞) has just 0, 1 more roots.
        q = Tλ(p, k)      ## shift p left by k so this is p on [k, ∞)

        no_left = descartes_bound(q)
        if no_left == 0
            break
        elseif no_left == 1
            _iszero(polyval(p, k)) ? push!(rts, k) : push!(rts, getIsolatedRoot(p, k, M+1))
            break
        end

        ## more to do, [k, ∞) might have 2 or more roots...
        st = find_in_01(q)
        for i in st.Exact
            push!(rts, k + i)
            p, mult = multiplicity(p, k+i)
        end
        for i in st.Isol
            push!(rts, k + getIsolatedRoot(q,i))
        end
    end
    rts
end
neg_roots(p::Poly) = pos_roots(polyval(p, -variable(p)))

## use Rational{BigInt} to find real roots
function real_roots{T <: QQ}(p::Poly{T})
    rts = Any[]
    p = convert(Poly{Rational{BigInt}}, p)
    p = divrem(p, bgcd(p, Polynomials.polyder(p)))[1] # square free p/gcd(p, p')
    p = p*(1//p[degree(p)])
    
    p,k = multiplicity(p, 0); k > 0 && push!(rts, 0) # 0
    append!(rts, pos_roots(p))                       # positive roots
    append!(rts, map(-,neg_roots(p)))                # negative roots
    rts
end


### Some functions to find all rational roots of a poly in Z[x], Q[x]
## return set of divisors of |n|, e.g. 12 -> [1,2,3,4,6,12]
function divisors(n::Integer)
    n == 0 && return([1])
    d = factor(abs(n))
    out = [1]
    for (k,v) in d
        out = unique(out * [k^i for i in 0:v]')
    end
    sort(out)
end

# find rational roots of p ∈ Z[x] 
function rational_roots{T <: Integer}(p::Poly{T})
    rts = Rational{T}[]
    if p[0] == 0
        push!(rts, 0)
        p, k = multiplicity(p, 0)
    else
        ds = divisors(p[0])
    end
    ns = divisors(p[degree(p)])
    poss = unique([d//n for d in ds, n in ns])
    for rt in poss
        for s in [-1,1]
            prt = s * rt
            polyval(p, prt) == 0 && push!(rts, prt)
        end
    end
    rts
end

# find rational roots of p ∈ Q[x] by using lcm(p)*p ∈ Z[x]
function rational_roots{T <: Integer}(p::Poly{Rational{T}})
    q = lcm(map(r -> r.den, p.a)) * p
    q = convert(Poly{BigInt}, q)
    rational_roots(q)
end




## factor over the rationals
## return (d,q) where d a dictionary of roots, and q divides p but has no rational factors.
function Base.factor{T <: QQ}(p::Poly{T})
    rts = rational_roots(p)
    d = Dict()
    for r in rts
        p,k = multiplicity(p,r)
        d[r] = k
    end
    d,p
end
