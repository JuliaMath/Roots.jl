using Base.Test
import Roots.fzero

# test problems from Table 1 of paper referenced in fzero.jl
# 1
@test_approx_eq fzero(x -> sin(x) - x/2, pi/2, pi) 1.89549426703398094714
@test_approx_eq Roots.a42a(x -> sin(x) - x/2, pi/2, pi) 1.89549426703398094714

# 2
function test2(x)
    -2*sum([(2i - 5)^2/(x - i^2)^3 for i = 1:20])
end
@test_approx_eq fzero(test2, 1 + 1e-9, 4 - 1e-9) 3.0229153472730568
@test_approx_eq_eps fzero(test2, 100 + 1e-9, 121 - 1e-9) 110.02653274766949 1e-9

# 3
function test3(x, a, b)
    a*x*exp(b*x)
end
@test_approx_eq fzero(x -> test3(x, -40, -1), -9, 31) 0.0
@test_approx_eq fzero(x -> test3(x, -100, -2), -9, 31) 0.0
@test_approx_eq fzero(x -> test3(x, -200, -3), -9, 31) 0.0

# 4
function test4(x, n, a)
    x^n - a
end
@test_approx_eq fzero(x -> test4(x, 4, 0.2), 0, 5) 0.668740304976422
@test_approx_eq fzero(x -> test4(x, 12, 1), 0, 5) 1.0
@test_approx_eq fzero(x -> test4(x, 8, 1), -0.95, 4.05) 1.0
@test_approx_eq fzero(x -> test4(x, 12, 1), -0.95, 4.05) 1.0

# 5
@test_approx_eq fzero(x -> sin(x) - 0.5, 0, 1.5) 0.5235987755982989

# 6
function test6(x, n)
    2x*exp(-n) - 2*exp(-n*x) + 1
end
@test_approx_eq fzero(x -> test6(x, 1), 0, 1) 0.42247770964123665
@test_approx_eq fzero(x -> test6(x, 20), 0, 1) 0.03465735902085385
@test_approx_eq fzero(x -> test6(x, 100), 0, 1) 0.006931471805599453

# 7
function test7(x, n)
    (1 + (1 - n)^2)*x - (1 - n*x)^2
end
@test_approx_eq fzero(x -> test7(x, 5), 0, 1) 0.0384025518406219
@test_approx_eq fzero(x -> test7(x, 10), 0, 1) 0.0099000099980005
@test_approx_eq fzero(x -> test7(x, 20), 0, 1) 0.0024937500390620117

# 8
function test8(x, n)
    x^2 - (1 - x)^n
end
@test_approx_eq fzero(x -> test8(x, 2), 0, 1) 0.5
@test_approx_eq fzero(x -> test8(x, 5), 0, 1) 0.345954815848242
@test_approx_eq fzero(x -> test8(x, 10), 0, 1) 0.24512233375330725
@test_approx_eq fzero(x -> test8(x, 20), 0, 1) 0.16492095727644096

# 9
function test9(x, n)
    (1 + (1 - n)^4)*x - (1 - n*x)^4
end
@test_approx_eq fzero(x -> test9(x, 1), 0, 1) 0.2755080409994844
@test_approx_eq fzero(x -> test9(x, 2), 0, 1) 0.1377540204997422
@test_approx_eq fzero(x -> test9(x, 20), 0, 1) 7.668595122185337e-6

# 10
function test10(x, n)
    exp(-n*x)*(x - 1) + x^n
end
@test_approx_eq fzero(x -> test10(x, 1), 0, 1) 0.401058137541547
@test_approx_eq fzero(x -> test10(x, 5), 0, 1) 0.5161535187579336
@test_approx_eq fzero(x -> test10(x, 20), 0, 1) 0.5527046666784878

# 11
function test11(x, n)
    (n*x - 1)/((n - 1)*x)
end
@test_approx_eq fzero(x -> test11(x, 2), 0.01, 1) 0.5
@test_approx_eq fzero(x -> test11(x, 15), 0.01, 1) 0.0666666666666666
@test_approx_eq fzero(x -> test11(x, 20), 0.01, 1) 0.05

# 12
function test12(x, n)
    x^(1/n) - n^(1/n)
end
@test_approx_eq fzero(x -> test12(x, 2), 1, 100) 2
@test_approx_eq fzero(x -> test12(x, 6), 1, 100) 6
@test_approx_eq fzero(x -> test12(x, 33), 1, 100) 33

# 13
function test13(x)
    if x == 0
        return 0
    else
        return x/exp(x^2)
    end
end
@test_approx_eq fzero(test13, -1, 4) 0

# 14
function test14(x, n)
    if x >= 0
        return n/20*(x/1.5 + sin(x) - 1)
    else
        return -n/20
    end
end
@test_approx_eq fzero(x -> test14(x, 1), -1e4, pi/2) 0.6238065189616123
@test_approx_eq fzero(x -> test14(x, 40), -1e4, pi/2) 0.6238065189616123

# 15
function test15(x, n)
    if x > 2e-3/(1+ n)
        return e - 1.859
    elseif x < 0
        return -0.859
    else
        return exp((n + 1)*x/2*1e3) - 1.859
    end
end
@test_approx_eq fzero(x -> test15(x, 20), -1e4, 1e4) 0.00005905130559421971
@test_approx_eq fzero(x -> test15(x, 40), -1e4, 1e4) 0.000030245790670210097
@test_approx_eq fzero(x -> test15(x, 100), -1e4, 1e4) 0.000012277994232461523
@test_approx_eq fzero(x -> test15(x, 1000), -1e4, 1e4) 1.2388385788997142e-6
