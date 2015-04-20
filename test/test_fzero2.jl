## Some more extreme tests from http://ir.igsnrr.ac.cn/bitstream/311030/8840/1/%E4%BE%AF%E9%BA%9F%E7%A7%91(SCI)2.pdf
## All the methods do poorly here due to the multiplicities of the functions
## even though there should be convergnce of either |xn+1 - xn| or |f(xn)|
using Roots
using Base.Test

tests = [
         (x -> (x - sqrt(5))^4 / ((x-1)^2 + 2),    4.3,  2.236067977499790)
         (x -> (sin(x)^2 - 2x + 1)^5,              4.5,  0.71483582544138924)
         (x -> (8x*exp(-x^2) -2x - 3)^8,          -2.0, -1.7903531791589544)
         (x -> (2x*cos(x) + x^2 - 3)^10/(x^2 + 1), 3.0,  2.9806452794385368)
         (x -> (exp(-x^2 + x + 3) - x + 2)^9,      2.5,  2.4905398276083051)
         (x -> (exp(-x) + 2sin(x))^4,              3.0,  3.1627488709263654)
         (x -> (log(x^2 + 3x + 5) - 2x + 7)^8,     5.5,  5.4690123359101421)
         (x -> (sqrt(x^2 + 2x + 5) - 
                           2sin(x) - x^2 + 3)^5,   2.3, 2.3319676558839640)
         (x -> (x-2)^4 / ( (x-1)^2 + 1),           1.8,  2.0000000000000000)
         (x -> abs(x - 2.5)^(15/4) * exp(x),       2.2,  2.5)
         (x -> (sqrt(x) - 1/x - 1)^7,              2.0,  2.147899035704787)
         (x -> (log(x) + sqrt(x) - 5)^3,           8.0,  8.309432694231572)
         (x -> (sin(x)*cos(x) - x^3 + 1)^9,        1.0,  1.117078770687451)
         (x -> ((x-3)*exp(x))^5,                   2.5,  3.0000000000000000)
         (x -> (log(x) + sqrt(x^4 + 1) -2)^7,      1.0,  1.222813963628973)
         ]


using Roots
ord(x) = string(floor(Integer,log10(x)))
orders = [0, 1, 2, 5, 8, 16]
out = Any[]
for (ctr, (f1, x0, xstar)) in enumerate(tests)
    m = String[]
    for order in orders
        try 
            a = fzero(f1, x0, order=order, maxeval=100)
            tm = begin tic(); fzero(f1, x0, order=order, maxeval=100); toc() end
            push!(m, ord(tm))
            push!(m, ord(abs(a-xstar)))
        catch
            push!(m, "***")
            push!(m, "***")
        end
    end
    try
        newton(f1, x0)
        tm =  begin tic(); a = newton(f1, x0); toc() end
        push!(m, ord(tm))
        push!(m, ord(abs(a-xstar)))
    catch e
        push!(m, "***")
        push!(m, "***")
    end
    push!(out, m)
end

## values in pairs time/accuracy or ****** if failure
## sensitive to maxeval value
 for i in out
     for k in 1:length(i)/2
         @printf "%3s%3s " i[2*(k-1)+1] i[2*k]
     end
         @printf "\n"
 end
