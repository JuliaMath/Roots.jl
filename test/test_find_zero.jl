## tests of find_zero interface
using Roots
using Compat.Test
import SpecialFunctions.erf


meths = [Order0(), Order1(), Order2(), Order5(), Order8(), Order16()]

## basics
fns = [(x -> x^2 - exp(x) - 3x + 2,   0.257530285439860, [-1, 0.5]),
       (x -> x^3 + 4x^2 -10,          1.365230013414097, [0.5]), # not -1.0, which causes errors in most
       (x -> exp(x)*sin(x) - 2x - 5, -2.523245230732555, [-1.0, 0.0]),
       (x -> log(x) + sqrt(x) - 5,    8.309432694231571, [15, 20]),
       (x -> sqrt(x) - 1/x - 3,       9.633595562832696, [20, 25]),

       (x -> exp(-x) - cos(x), 1.292695719373398, [-1.0, .75, 1.5]),  # not 0.6 which has issues with the bump
       (x -> cos(x)^2 - x/5, 1.085982678007472, [0.5, 1.0]),
       (x -> x^10 - x^3 - x - 1, -0.674177935277052, [0.25, 0.5]),
       (x -> sin(x) - x + 1, 1.934563210752024, [2.0, 20.0]),
       (x -> exp(-x^2 + x + 2) - cos(x) + x^3 + 1, -1.0899423334391694, [-0.9, -0.6]),

       (x -> asin(x^2 -1) - x/2 + 1, 0.594810968398369, [0.3, 0.9]),
       (x -> tanh(x) - tan(x), 7.068582745628732, [5.5, 7.1])   # not 5.0, as some methods go off to vertical asymptote
]


for (i, (f, xstar, xs)) in enumerate(fns)
    for m in meths
        for x0_ in xs
            out = try
                xn = find_zero(f, x0_, m)
                @test abs(xn - xstar) < 1e-14 || abs(f(xn)) < 1e-13
                "."
            catch err
                "*"
            end
            out == "*" && println("Fn $i, method $m, x0=$x0_ failed")
        end
    end
end


## Some more extreme tests from http://ir.igsnrr.ac.cn/bitstream/311030/8840/1/%E4%BE%AF%E9%BA%9F%E7%A7%91(SCI)2.pdf
## All the methods do poorly here due to the multiplicities of the functions
## even though there should be convergnce of either |xn+1 - xn| or |f(xn)|
multiplicity_tests = [
         (x -> (x - sqrt(5))^4 / ((x-1)^2 + 2),    4.3,  2.236067977499790)
         (x -> (sin(x)^2 - 2x + 1)^5,              4.5,  0.71483582544138924)
         (x -> (8x*exp(-x^2) -2x - 3)^8,          -2.0, -1.7903531791589544)
         (x -> (2x*cos(x) + x^2 - 3)^10/(x^2 + 1), 3.0,  2.9806452794385368)
         (x -> (exp(-x^2 + x + 3) - x + 2)^9,      2.5,  2.4905398276083051)

         (x -> (exp(-x) + 2sin(x))^4,              3.0,  3.1627488709263654)
         (x -> (log(x^2 + 3x + 5) - 2x + 7)^8,     5.5,  5.4690123359101421)
         (x -> (sqrt(x^2 + 2x + 5) - 2sin(x) - x^2 + 3)^5,   2.3, 2.3319676558839640)
         (x -> (x-2)^4 / ( (x-1)^2 + 1),           1.8,  2.0000000000000000)
         (x -> abs(x - 2.5)^(15/4) * exp(x),       2.2,  2.5)

         (x -> (sqrt(x) - 1/x - 1)^7,              2.0,  2.147899035704787)
         (x -> (log(x) + sqrt(x) - 5)^3,           8.0,  8.309432694231572)
         (x -> (sin(x)*cos(x) - x^3 + 1)^9,        1.0,  1.117078770687451)
         (x -> ((x-3)*exp(x))^5,                   2.5,  3.0000000000000000)
         (x -> (log(x) + sqrt(x^4 + 1) -2)^7,      1.0,  1.222813963628973)
         ]

for (i, (fn_, x0_, xstar)) in enumerate(multiplicity_tests)
    for m in meths
        #println("$i: $m")
        @test  abs(find_zero(fn_, x0_, m, maxevals=100) - xstar) < 1e-1 # wow, not too ambitious here, 9th powers...
    end
end


## issues with starting near a maxima. Some bounce out of it, but
## one would expect all to have issues
fn, xstar = x -> x^3 + 4x^2 -10,  1.365230013414097
for m in [Order0(), Order1()]
    @test_throws Roots.ConvergenceFailed find_zero(fn, -1.0, m)
end
for m in [Order2(), Order5(), Order8(), Order16()]
    @test find_zero(fn, -1.0, m) ≈ xstar
end

## non-simple root
fn, xstar, x0 = x -> cos(x) - 1, 0.0, 0.1
for m in meths
    xn = find_zero(fn, x0, m)
    @test abs(fn(xn)) <= 1e-10
end

## issue with large steps
fn, x0 = x -> x^20 - 1, 0.5
for m in meths[2:end] # not 0, as it uses bracket
    @test_throws Roots.ConvergenceFailed find_zero(fn, x0, m)
end


## different types of initial values
for m in meths
    @test find_zero(sin, 3, m) ≈ pi
    @test find_zero(sin, 3.0, m) ≈ pi
    @test find_zero(sin, big(3), m) ≈ pi
    @test find_zero(sin, big(3.0), m) ≈ pi
end


## Methods guarded with a bracket
# fn, xstar, x0 = (x -> sin(x) - x - 1,  -1.9345632107520243, 2)
# @test_throws Roots.ConvergenceFailed find_zero(fn, x0, Order2())
# for m in meths
#     @test find_zero(fn, x0, m, bracket=[-2,3]) ≈ xstar
# end


# fn, xstar, x0 = (x -> x * exp( - x ), 0, 1.0)
# @test find_zero(fn, x0, Order0(), bracket=[-1,2]) ≈ xstar
# @test find_zero(fn, 7.0, Order0(), bracket=[-1,2]) ≈ xstar  # out of bracket

## bisection methods
@test find_zero(x -> cos(x) - x, [0, pi], Bisection()) ≈ 0.7390851332151607
@test (find_zero(x -> sin(x), [big(3), 4], Bisection()) |> sin) < 1e-70
@test find_zero(x -> sin(x), [4,3], Bisection()) ≈ pi
@test find_zero(x -> sin(x), [big(4),3], Bisection()) ≈ pi
@test find_zero(x -> Inf*sign(x-1.0), (-Inf, Inf)) ≈ 1
#@test_throws ArgumentError  find_zero(x -> sin(x), 3.0, Bisection())

@test find_zero(x -> cos(x) - x, [0, pi], FalsePosition()) ≈ 0.7390851332151607
@test (find_zero(x -> sin(x), [big(3), 4], FalsePosition(), maxevals=1000) |> sin) < 1e-70
@test find_zero(x -> sin(x), [4,3], FalsePosition()) ≈ pi

galadino_probs = [(x -> x^3 - 1, [.5, 1.5]),
                  (x -> x^2 * (x^2/3 + sqrt(2) * sin(x)) - sqrt(3)/18, [.1, 1]),
                  (x -> 11x^11 - 1, [0.1, 1]),
                  (x ->  x^3 + 1, [-1.8, 0]),
                  (x ->  x^3 - 2x - 5, [2.0, 3]),
                  
                  ((x,n=5)  -> 2x * exp(-n) + 1 - 2exp(-n*x) , [0,1]),
                  ((x,n=10) -> 2x * exp(-n) + 1 - 2exp(-n*x) , [0,1]),
                  ((x,n=20) -> 2x * exp(-n) + 1 - 2exp(-n*x) , [0,1]),

                  ((x,n=5)  -> (1 + (1-n)^2) * x^2 - (1 - n*x)^2 , [0,1]),
                  ((x,n=10) -> (1 + (1-n)^2) * x^2 - (1 - n*x)^2 , [0,1]),
                  ((x,n=20) -> (1 + (1-n)^2) * x^2 - (1 - n*x)^2 , [0,1]),

                  ((x,n=5)  -> x^2 - (1-x)^n , [0,1]),
                  ((x,n=10) -> x^2 - (1-x)^n , [0,1]),
                  ((x,n=20) -> x^2 - (1-x)^n , [0,1]),

                  ((x,n=5)  -> (1 + (1-n)^4)*x - (1 - n*x)^4 , [0,1]),
                  ((x,n=10) -> (1 + (1-n)^4)*x - (1 - n*x)^4 , [0,1]),
                  ((x,n=20) -> (1 + (1-n)^4)*x - (1 - n*x)^4 , [0,1]),

                  ((x,n=5)  -> exp(-n*x)*(x-1) + x^n , [0,1]),
                  ((x,n=10) -> exp(-n*x)*(x-1) + x^n , [0,1]),
                  ((x,n=20) -> exp(-n*x)*(x-1) + x^n , [0,1]),

                  ((x,n=5)  -> x^2 + sin(x/n) - 1/4 , [0,1]),
                  ((x,n=10) -> x^2 + sin(x/n) - 1/4 , [0,1]),
                  ((x,n=10) -> x^2 + sin(x/n) - 1/4 , [0,1])
                  ]
        

for (fn_, ab) in galadino_probs
    for m in [Roots.A42(), Bisection(), (FalsePosition(i) for i in 1:12)...]
        global x0 = find_zero(fn_, ab, m, maxevals=120)
        @test abs(fn_(x0)) <= 1e-14
    end
end

## Can use (a,b) or [a,b] for bisection
fn = x -> x^5 - x - 1
for m in [Roots.A42(), Bisection(), (FalsePosition(i) for i in 1:12)...]
    global x0 = find_zero(fn, (1,2.0), m)
    @test abs(fn(x0)) <= 1e-14
end


## defaults for method argument
@test find_zero(x -> cbrt(x), 1) ≈ 0.0 # order0()
@test find_zero(sin, [3,4]) ≈ π   # Bisection() 


## Callable objects
struct _SampleCallableObject end
_SampleCallableObject(x) = x^5 - x - 1
for m in meths
    @test find_zero(_SampleCallableObject, 1.0, m) ≈ 1.1673039782614187
end

### a wrapper to count function calls, say
mutable struct Cnt
    cnt::Int
    f
    Cnt(f) = new(0, f)
end
(f::Cnt)(x) = (f.cnt += 1; f.f(x))

g = Cnt(x -> x^5 - x - 1)
for m in meths
    @test find_zero(g, 1.0, m) ≈ 1.1673039782614187
end


## test tolerance arguments
fn, xstar = x -> sin(x) - x + 1, 1.9345632107520243
@test find_zero(fn, 20.0, Order2())  ≈ xstar   # needs 16 iterations, 33 fn evaluations
@test abs(fn(find_zero(fn, 20.0, Order2(), atol=1e-2)) - xstar) > 1e-12
@test abs(fn(find_zero(fn, 20.0, Order2(), rtol=1e-2)) - xstar) > 1e-12
@test_throws Roots.ConvergenceFailed find_zero(fn, 20.0, Order2(), maxevals=5) 
@test_throws Roots.ConvergenceFailed find_zero(fn, 20.0, Order2(), maxfnevals=10) 


## test robustness of Order0
## tests from: http://people.sc.fsu.edu/~jburkardt/cpp_src/test_zero/test_zero.html
function _newton_baffler(x) 
    if ( x - 0.0 ) < -0.25 
        0.75 * ( x - 0 ) - 0.3125 
    elseif  ( x - 0 ) < 0.25 
        2.0 * ( x - 0 ) 
    else
        0.75 * ( x - 0 ) + 0.3125
    end
end
pathological = [
                (x ->sin( x ) - x / 2, .1),
                (x-> 2 * x - exp( - x ), 1),
                (x -> x * exp( - x ), 0.99),
                (x -> exp( x ) - 1 / ( 10 * x )^2, .1),
                (x -> ( x + 3 ) * ( x - 1 )^2, 1),
                
                (x -> exp( x ) - 2 - 1 / ( 10 * x )^2 + 2 / ( 100 * x )^3, 1),
                (x -> cbrt(x), 1),
                (x -> cos( x ) - x, 1),
                (x-> _newton_baffler(x), 8),
                (x -> 20.0 * x / ( 100.0 * x^2 + 1.0), 0.095), 
                
                (x ->  ( 4.0 + x^2) * ( 2.0 + x ) * ( 2.0 - x )  / ( 16.0 * x * x * x * x + 0.00001 ), 1),
                (x -> (x == 1.0) ? float(0) : sign(x-1.0) * exp(log(1e4) + log(abs(x - 1.0)) - 1.0/(x-1.0)^2), 1),
                (x -> 0.00000000001 * (x - 100.0), 1),
                (x -> 1.0 / ( ( x - 0.3 ) * ( x - 0.3 ) + 0.01 ) + 1.0 / ( ( x - 0.9 ) * ( x - 0.9 ) + 0.04 ) + 2.0 * x - 5.2, -1),
                (x -> ( 1 - 6x^2) * cbrt(x) * exp(-x^2) / (3*x), -0.25), 
                
                (x -> ( pi * ( x - 5.0 ) / 180.0 ) - 0.8 * sin( pi * x / 180.0 ), 1),
                (x -> x^3 - 2*x - 5, 2),
(x -> 1e6 * (x^7 -7x^6 +21x^5 -35x^4 +35x^3-21x^2+7x-1),  0.990),
(x -> cos(100*x)-4*erf(30*x-10), 0.0) 
                ]


for (fn_, x0_) in pathological
    find_zero(fn_, x0_, Order0())
end

## issue tests: put in tests to ensure closed issues don't reappear.

## issue #94; tolerances not matching documentation
test_94 = function(;kwargs...)
    g, T = 1.62850, 14.60000
    α, t1, tf = 0.00347, 40.91375, 131.86573
    y, ya, yf = 0.0, 9000.0, 10000.0
    vy = sqrt(2g*(ya-y))
    θ0, θ1 = atan(α*tf), atan(α*(tf-t1))
    I_sintan(x) = tan(x)/2cos(x) - atanh(tan(x/2))
    I_sintan(x, y) = I_sintan(y) - I_sintan(x)
    function lhs(θ)
        tRem = (vy - T/α*(sec(θ1) - sec(θ))) / g
        val = -yf + y + vy*tRem - 0.5g*tRem^2 - T/α^2*I_sintan(θ, θ1)
        val
    end

    meth = Roots.FalsePosition()
    f, x0 = lhs, [atan(α*tf), atan(α*(tf-t1))]
    F = Roots.DerivativeFree(lhs)
    state = Roots.init_state(meth, F, x0)
    options = Roots.init_options(meth, state)
    find_zero(meth, F, options, state)

    @test state.steps <= 15
end
test_94()

@testset "init_state!" begin
    ## add init_state! method
    g1(x) = x^5 - x - 1
    x0_, xstar_ = (1.0, 2.0), 1.1673039782614187
    M = Roots.A42()
    state = Roots.init_state(M, g1, x0_)
    options = Roots.init_options(M, state)
    for M in (Roots.A42(), Roots.Bisection(), Roots.FalsePosition())
        Roots.init_state!(state, M, g1, x0_)
        @test find_zero(M, g1, options, state) ≈ xstar_
    end


    ## hybrid
    g1(x) = exp(x) - x^4
    x0_, xstar_ = (5.0, 20.0), 8.613169456441398
    M = Roots.Bisection()
    state = Roots.init_state(M, g1, x0_)
    options = Roots.init_options(M, state, xatol=1/2)
    find_zero(M, g1, options, state)
    M = Roots.Order1()
    Roots.init_state!(state, M, g1, state.m[1])
    Roots.init_options!(options, M)
    @test find_zero(M, g1, options, state) ≈ xstar_


    # test conversion between states/options
    Ms = vcat([Roots.FalsePosition(i) for i in 1:12], Roots.A42(), Roots.AlefeldPotraShi(), Roots.Brent(), Roots.Bisection())
    Ns = [Roots.Order1(), Roots.Order2(), Roots.Order5(), Roots.Order8(), Roots.Order16()]

    g1(x) = exp(x) - x^4
    x0_, xstar_ = (5.0, 20.0), 8.613169456441398
    M = Ms[1]
    state = Roots.init_state(M, g1, x0_)
    options = Roots.init_options(M, state)

    for M in Ms
        @test Roots.init_state!(state, M, g1, x0_) == nothing
        @test Roots.init_options!(options, M) == nothing
        @test Roots.update_state(M, g1, state, options) == nothing
        @test Roots.assess_convergence(M, state, options) == false # none in just 1 step
    end

    x0_ = 8.0
    for M in Ns
        @test Roots.init_state!(state, M, g1, x0_) == nothing
        @test Roots.init_options!(options, M) == nothing
        @test Roots.update_state(M, g1, state, options)  == nothing
        @test Roots.assess_convergence(M, state, options) == false  # none in 1 step
    end  
    
end
    
    
    
    

