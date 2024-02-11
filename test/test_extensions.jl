#=
using SymPy
@testset "SymPy" begin
    SymPy.@syms x

    @test find_zero(cos(x) ~ 1/2, (0, pi/2)) ≈ find_zero(x -> cos(x) - 1/2, (0, pi/2))
    @test find_zero(1/2 ~ cos(x), (0, pi/2)) ≈ find_zero(x -> 1/2 - cos(x), (0, pi/2))
    @test find_zero(cos(x) ~ x/2, (0, pi/2)) ≈ find_zero(x -> cos(x) - x/2, (0, pi/2))

    @test find_zeros(cos(x) ~ 1/2, (0, pi/2)) ≈ find_zeros(x -> cos(x) - 1/2, (0, pi/2))
    @test find_zeros(1/2 ~ cos(x), (0, pi/2)) ≈ find_zeros(x -> 1/2 - cos(x), (0, pi/2))
    @test find_zeros(cos(x) ~ x/2, (0, pi/2)) ≈ find_zeros(x -> cos(x) - x/2, (0, pi/2))

    @test fzero(cos(x) ~ 1/2, 0, pi/2) ≈ fzero(x -> cos(x) - 1/2, 0, pi/2)
    @test fzero(1/2 ~ cos(x), 0, pi/2) ≈ fzero(x -> 1/2 - cos(x), 0, pi/2)
    @test fzero(cos(x) ~ x/2, 0, pi/2) ≈ fzero(x -> cos(x) - x/2, 0, pi/2)

    @test fzeros(cos(x) ~ 1/2, 0, pi/2) ≈ fzeros(x -> cos(x) - 1/2, 0, pi/2)
    @test fzeros(1/2 ~ cos(x), 0, pi/2) ≈ fzeros(x -> 1/2 - cos(x), 0, pi/2)
    @test fzeros(cos(x) ~ x/2, 0, pi/2) ≈ fzeros(x -> cos(x) - x/2, 0, pi/2)

end
=#

#=
using SymPyPythonCall
@testset "SymPythonCall" begin
    SymPyPythonCall.@syms x

    @test find_zero(cos(x) ~ 1/2, (0, pi/2)) ≈ find_zero(x -> cos(x) - 1/2, (0, pi/2))
    @test find_zero(1/2 ~ cos(x), (0, pi/2)) ≈ find_zero(x -> 1/2 - cos(x), (0, pi/2))
    @test find_zero(cos(x) ~ x/2, (0, pi/2)) ≈ find_zero(x -> cos(x) - x/2, (0, pi/2))

    @test find_zeros(cos(x) ~ 1/2, (0, pi/2)) ≈ find_zeros(x -> cos(x) - 1/2, (0, pi/2))
    @test find_zeros(1/2 ~ cos(x), (0, pi/2)) ≈ find_zeros(x -> 1/2 - cos(x), (0, pi/2))
    @test find_zeros(cos(x) ~ x/2, (0, pi/2)) ≈ find_zeros(x -> cos(x) - x/2, (0, pi/2))

    @test fzero(cos(x) ~ 1/2, 0, pi/2) ≈ fzero(x -> cos(x) - 1/2, 0, pi/2)
    @test fzero(1/2 ~ cos(x), 0, pi/2) ≈ fzero(x -> 1/2 - cos(x), 0, pi/2)
    @test fzero(cos(x) ~ x/2, 0, pi/2) ≈ fzero(x -> cos(x) - x/2, 0, pi/2)

    @test fzeros(cos(x) ~ 1/2, 0, pi/2) ≈ fzeros(x -> cos(x) - 1/2, 0, pi/2)
    @test fzeros(1/2 ~ cos(x), 0, pi/2) ≈ fzeros(x -> 1/2 - cos(x), 0, pi/2)
    @test fzeros(cos(x) ~ x/2, 0, pi/2) ≈ fzeros(x -> cos(x) - x/2, 0, pi/2)

end
=#

using ForwardDiff
@testset "ForwardDiff" begin
    f = (x, p) -> x^2 - p
    Z = ZeroProblem(f, (0, 1000))
    F = p -> solve(Z, Roots.Bisection(), p)
    for p in (3, 5, 7, 11)
        @test F(p) ≈ sqrt(p)
        @test ForwardDiff.derivative(F, p) ≈ 1 / (2sqrt(p))
    end

    # Hessian is *fixed* for F(p) = find_zero(f, x₀, M, p)
    f = (x, p) -> x^2 - sum(p .^ 2)
    Z = ZeroProblem(f, (0, 1000))
    F = p -> solve(Z, Roots.Bisection(), p)
    Z = ZeroProblem(f, (0, 1000))
    F = p -> solve(Z, Roots.Bisection(), p)
    for p in ([1, 2], [1, 3], [1, 4])
        @test F(p) ≈ sqrt(sum(p .^ 2))
        a, b = p
        n = sqrt(a^2 + b^2)^3
        @test ForwardDiff.hessian(F, p) ≈ [b^2 -a*b; -a*b a^2] / n
    end
end

#=
using IntervalRootFinding
@testset "IntervalRootFinding" begin
    f(x) = sin(x + sin(x + sin(x)))
    @test find_zeros(f, (-5, 5)) ≈ [-pi, 0, pi]
    out = find_zeros(f, -5..5, Roots.Newton())
    @test sort(out.zeros) ≈ sort([-pi,0,pi])
    @test isempty(out.unknown)
end
=#
