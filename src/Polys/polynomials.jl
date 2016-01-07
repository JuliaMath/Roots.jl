## Extensions to Polynomials.jl

" Return the variable `x` of a polynomial as a polynomial."
variable(p::Poly) = poly(zeros(eltype(p),1), p.var)

"Create a monic polynomial from `p`"
monic(p::Poly) = Poly(p.a/p[degree(p)], p.var)

## Reverse coefficients
rcoeffs(p::Poly) = reverse(coeffs(p))


function Base.convert{T<:Integer}(::Type{Poly{T}}, p::Poly{Rational{T}})
    l = reduce(lcm, [x.den for x in p.a])
    q = l * p
    Poly(T[x for x in q.a], p.var)
end

## convert Poly <--> function
Base.convert(::Type{Function}, p::Poly) = x -> Polynomials.polyval(p,x)
## convert a function to a polynomial with error if conversion is not possible
## This needs a hack, as Polynomials.jl defines `/` to return a `div` and not an error for p(x)/q(x)
immutable PolyTest                                                               
   x                                                                             
end                                                                              
import Base: +,-,*,/,^                                                           
+(a::PolyTest, b::PolyTest) = PolyTest(a.x + b.x)                                
+{T<:Number}(a::T, b::PolyTest) = PolyTest(a + b.x)                              
+{T<:Number}(a::PolyTest,b::T) = PolyTest(a.x + b)                              
-(a::PolyTest,b::PolyTest) = PolyTest(a.x - b.x)                                
-{T<:Number}(a::T, b::PolyTest) = PolyTest(a - b.x)                              
-{T<:Number}(a::PolyTest,b::T) = PolyTest(a.x - b)                              
-(a::PolyTest) = PolyTest(-a.x)                                                 
*(a::PolyTest, b::PolyTest) = PolyTest(a.x * b.x)                               
*(a::Bool, b::PolyTest) = PolyTest(a * b.x)                                     
*{T<:Number}(a::T, b::PolyTest) = PolyTest(a * b.x)                             
*{T<:Number}(b::PolyTest, a::T) = PolyTest(a * b.x)                             
/{T<:Number}(b::PolyTest, a::T) = PolyTest(b.x / a)                             
^{T<:Integer}(b::PolyTest, a::T) = PolyTest(b.x ^ a)                            

typealias QQR @compat Union{Int, BigInt, Rational{Int}, Rational{BigInt}, Float64}
function Base.convert{T<:QQR}(::Type{Poly{T}}, f::Function)
    try
        f(PolyTest(0))                  # error if not a poly
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


