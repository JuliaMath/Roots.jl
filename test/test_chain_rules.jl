using Roots
using Zygote
using Test

# issue #325 add frule, rrule

@testset "Test rrule" begin

    # single function
    f(x, p) = log(x) - p
    F(p) = find_zero(f, 1, Order1(), p)
    @test first(Zygote.gradient(F, 1)) ≈ exp(1)

    g(x, p) = x^2 - p[1]*x - p[2]
    G(p) = find_zero(g, 1, Order1(), p)
    @test first(Zygote.gradient(G, [0,4])) ≈ [1/2, 1/4]

    fx(x,p) = 1/x
    F2(p) = find_zero((f,fx), 1, Roots.Newton(), p)
    @test first(Zygote.gradient(F2, 1)) ≈ exp(1)

    gp(x, p) = 2x - p[1]
    G2(p) = find_zero((g, gp), 1, Roots.Newton(), p)
    @test first(Zygote.gradient(G2, [0,4])) ≈ [1/2, 1/4]

end
