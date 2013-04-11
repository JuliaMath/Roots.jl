using Roots

root = fzero(sin, -0.5, 0.5)
println(sin(root))

f(x) = (x - 1)^3
root = fzero(f, 0.5, 2.0)
println(f(root))

f(x) = 2x*exp(-20) - 2*exp(-20x) + 1

roots = fzero(f, 0.0, 1.0; tolerance=1e-10, max_iter=100)
println(f(root))


f(x) = exp(x) - cos(x)
fp(x) = exp(x) + sin(x)

root = newton(f, fp, 3.0)
println(f(root))


f(x) = exp(x) - cos(x)
fp(x) = exp(x) + sin(x)
fpp(x) = exp(x) + cos(x)

root = halley(f, fp, fpp, 3.0)
println(f(root))

using Calculus

f(x) = exp(x) - cos(x)
root = newton(f, derivative(f), 3.0)
println(f(root))

