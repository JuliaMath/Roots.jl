##module AD

## Forward automatic differentiation code
## Basic use is for f:R -> R
##
## Unlike symbolic differentiation, this exports an operator D which
## maps a function to a function (not an operator that acts on expressions).
## f(x) = exp(-x)*sin(x)
## fp(x) = D(f)(x) ...
## fpp(x) = D(f, 2)(x) ...

importall Base

export D, D2

immutable Ad
    val
    der
end

ad(x::Real) = Ad(x, 1)
ad{T <: Real}(x::Array{T}) = Ad(x, ones(length(x)))
ad(x::Ad) = Ad(x, 1)
ad(x::Array{Ad}) = Ad(x, ones(length(x)))
 

## conversion
convert(::Type{Ad}, x::Real) = Ad(x, 0)
promote_rule{T <: Real}(::Type{Ad}, ::Type{T}) = Ad


## accessors
val(a::Ad) = a.val
val(a::Real) = a
der(a::Ad) = a.der
der(a::Real) = 0

## Need to extract
function ad_order(a::Ad)
    local b = a
    ctr = 1
    while isa(b.val, Ad)
        ctr += 1
        b = b.val
    end
    ctr
end

## compose functions with *
*(f::Function, g::Function) = x -> f(g(x))
getindex(x::Ad,i::Int) = x.val[i]

## compute f^(k)(a)
function Base.diff(a::Ad, k::Integer)
    
    if k == 0 return(a) end
    n = ad_order(a)
    if k > n
        error("Can only take upto $n derivatives for this value of a")
    end

    if k == n
      (der^k)(a)
    else
      (val^(n-k))((der^k)(a))
    end
end

    
## Operator for derivatives
## 
## @param f: R -> R
## @param k degree of derivative, defaults to 1
function D(f::Function, k::Integer)
    if k == 0
        x -> f(x)
    elseif k > 0
        x -> diff(f((ad^k)(x)), k)
    else
        stop("Can't handle negative k")
    end
end
D(fs::Array{Function}) = map(D, fs)


D(f::Function) = D(f, 1)
D2(f::Function) = D(f, 2)




## math ops
+(x::Ad, y::Ad)           = Ad(x.val + y.val, x.der + y.der)
+{T <: Real}(x::Ad, y::T) = +(promote(x,y)...)
+{T <: Real}(x::T, y::Ad) = +(promote(x,y)...)

.+(x::Array{Ad}, y::Array{Ad})           = Ad(x.val .+ y.val, x.der .+ y.der)
.+{T <: Real}(x::Array{Ad}, y::Array{T}) = .+(promote(x,y)...)
.+{T <: Real}(x::Array{T}, y::Array{Ad}) = .+(promote(x,y)...)


-(x::Ad) = Ad(-x.val, -x.der)


-(x::Ad, y::Ad)           = Ad(x.val - y.val, x.der - y.der)
-{T <: Real}(x::Ad, y::T) = -(promote(x,y)...)
-{T <: Real}(x::T, y::Ad) = -(promote(x,y)...)

.-(x::Array{Ad}, y::Array{Ad})           = Ad(x.val .- y.val, x.der .- y.der)
.-{T <: Real}(x::Array{Ad}, y::Array{T}) = .-(promote(x,y)...)
.-{T <: Real}(x::Array{T}, y::Array{Ad}) = .-(promote(x,y)...)


*(x::Ad, y::Ad)           = Ad(x.val * y.val, x.val * y.der + x.der * y.val)
*{T <: Real}(x::Ad, y::T) = *(promote(x,y)...)
*{T <: Real}(x::T, y::Ad) = *(promote(x,y)...)

.*(x::Array{Ad}, y::Array{Ad})           = Ad(x.val .* y.val, x.val .* y.der .+ x.der .* y.val)
.*{T <: Real}(x::Array{Ad}, y::Array{T}) = .*(promote(x,y)...)
.*{T <: Real}(x::Array{T}, y::Array{Ad}) = .*(promote(x,y)...)


/(x::Ad, y::Ad)           = Ad(x.val / y.val, (x.der * y.val - x.val * y.der)/(y.val^2))
/{T <: Real}(x::Ad, y::T) = /(promote(x,y)...)
/{T <: Real}(x::T, y::Ad) = /(promote(x,y)...)

./(x::Array{Ad}, y::Array{Ad})           = Ad(x.val ./ y.val, (x.der .* y.val .- x.val .* y.der)./(y.val.^2))
./{T <: Real}(x::Array{Ad}, y::Array{T}) = ./(promote(x,y)...)
./{T <: Real}(x::Array{T}, y::Array{Ad}) = ./(promote(x,y)...)


one(x::Ad) = Ad(1, 0)                        ## multiplicative identity
function ^(x::Ad, y::Integer) 
    if y == 0
        one(x)
    elseif y < 0 
        (1/x)^abs(y) 
    else 
        Ad(x.val ^ y, y * (x.val) ^ (y-1) * x.der) 
    end
end

^(x::Ad, y::Ad)   = exp(y * log(x))
^(x::Ad, y::Real) = Ad(x.val ^ y, y * (x.val) ^ (y-1) * x.der)
^(x::Real, y::Ad) = exp(y * log(x))

.^(x::Ad, y::Real) = Ad(x.val .^ y, y .* (x.val) .^ (y.-1) .* x.der)
.^(x::Ad, y::Ad)   = exp(y .* log(x))
.^(x::Real, y::Ad) = exp(y .* log(x))

## comparison
isless(a::Ad, y::Real) = isless(a.val, y)
isless(x::Real, b::Ad) = isless(x, b.val)
isless(a::Ad, b::Ad)   = isless(a.val, b.val)

for k in (:<, :<=, :>=, :>, :(==))
    @eval begin
        ($k)(x::Ad, y::Ad) = ($k)(x.val, y.val)
    end
end


for k in (:max, :min)
    @eval begin
        ($k)(x::Ad, y::Ad) = ($k)(x.val, y.val)
    end
end

## From the Calculus package
Calculus_derivative_rules = [
    ( :sqrt,        :(  xp / 2 / sqrt(x)                         ))
    ( :cbrt,        :(  xp / 3 / cbrt(x)^2                       ))
    ( :square,      :(  xp * 2 * x                               ))
    ( :log,         :(  xp / x                                   ))
    ( :log10,       :(  xp / x / log(10)                         ))
    ( :log2,        :(  xp / x / log(2)                          ))
    ( :log1p,       :(  xp / (x + 1)                             ))
    ( :exp,         :(  xp * exp(x)                              ))
    ( :exp2,        :(  xp * log(2) * exp2(x)                    ))
    ( :expm1,       :(  xp * exp(x)                              ))
    ( :sin,         :(  xp * cos(x)                              ))
    ( :cos,         :( -xp * sin(x)                              ))
    ( :tan,         :(  xp * (1 + tan(x)^2)                      ))
    ( :sec,         :(  xp * sec(x) * tan(x)                     ))
    ( :csc,         :( -xp * csc(x) * cot(x)                     ))
    ( :cot,         :( -xp * (1 + cot(x)^2)                      ))
    ( :sind,        :(  xp * cosd(x)                             ))
    ( :cosd,        :( -xp * sind(x)                             ))
    ( :tand,        :(  xp * (1 + tand(x)^2)                     ))
    ( :secd,        :(  xp * secd(x) * tand(x)                   ))
    ( :cscd,        :( -xp * cscd(x) * cotd(x)                   ))
    ( :cotd,        :( -xp * (1 + cotd(x)^2)                     ))
    ( :asin,        :(  xp / sqrt(1 - x^2)                       ))
    ( :acos,        :( -xp / sqrt(1 - x^2)                       ))
    ( :atan,        :(  xp / (1 + x^2)                           ))
    ( :asec,        :(  xp / abs(x) / sqrt(x^2 - 1)              ))
    ( :acsc,        :( -xp / abs(x) / sqrt(x^2 - 1)              ))
    ( :acot,        :( -xp / (1 + x^2)                           ))
    ( :asind,       :(  xp * 180 / pi / sqrt(1 - x^2)            ))
    ( :acosd,       :( -xp * 180 / pi / sqrt(1 - x^2)            ))
    ( :atand,       :(  xp * 180 / pi / (1 + x^2)                ))
    ( :asecd,       :(  xp * 180 / pi / abs(x) / sqrt(x^2 - 1)   ))
    ( :acscd,       :( -xp * 180 / pi / abs(x) / sqrt(x^2 - 1)   ))
    ( :acotd,       :( -xp * 180 / pi / (1 + x^2)                ))
    ( :sinh,        :(  xp * cosh(x)                             ))
    ( :cosh,        :(  xp * sinh(x)                             ))
    ( :tanh,        :(  xp * sech(x)^2                           ))
    ( :sech,        :( -xp * tanh(x) * sech(x)                   ))
    ( :csch,        :( -xp * coth(x) * csch(x)                   ))
    ( :coth,        :( -xp * csch(x)^2                           ))
    ( :asinh,       :(  xp / sqrt(x^2 + 1)                       ))
    ( :acosh,       :(  xp / sqrt(x^2 - 1)                       ))
    ( :atanh,       :(  xp / (1 - x^2)                           ))
    ( :asech,       :( -xp / x / sqrt(1 - x^2)                   ))
    ( :acsch,       :( -xp / abs(x) / sqrt(1 + x^2)              ))
    ( :acoth,       :(  xp / (1 - x^2)                           ))
    ( :erf,         :(  xp * 2 * exp(-square(x)) / sqrt(pi)      ))
    ( :erfc,        :( -xp * 2 * exp(-square(x)) / sqrt(pi)      ))
    ( :erfi,        :(  xp * 2 * exp(square(x)) / sqrt(pi)       ))
    ( :gamma,       :(  xp * digamma(x) * gamma(x)               ))
    ( :lgamma,      :(  xp * digamma(x)                          ))
    ( :airy,        :(  xp * airyprime(x)                        ))  # note: only covers the 1-arg version
    ( :airyprime,   :(  xp * airy(2, x)                          ))
    ( :airyai,      :(  xp * airyaiprime(x)                      ))
    ( :airybi,      :(  xp * airybiprime(x)                      ))
    ( :airyaiprime, :(  xp * x * airyai(x)                       ))
    ( :airybiprime, :(  xp * x * airybi(x)                       ))
    ( :besselj0,    :( -xp * besselj1(x)                         ))
    ( :besselj1,    :(  xp * (besselj0(x) - besselj(2, x)) / 2   ))
    ( :bessely0,    :( -xp * bessely1(x)                         ))
    ( :bessely1,    :(  xp * (bessely0(x) - bessely(2, x)) / 2   ))
    ## ( :erfcx,   :(  xp * (2 * x * erfcx(x) - 2 / sqrt(pi))   ))  # uncertain
    ## ( :dawson,  :(  xp * (1 - 2x * dawson(x))                ))  # uncertain

                             ]


for (e, ex) in Calculus_derivative_rules

    @eval begin
        ($e)(x::Ad) = Ad($(e)(x.val), 
        begin
            tmp = @eval((x, xp) -> ($ex))
            tmp(x.val, x.der)
        end)
    end
end
