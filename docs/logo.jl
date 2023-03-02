using Plots, Roots, FileIO
julia_colors = [:royalblue, :brown3, :forestgreen,  :mediumorchid3]

f(x) =  x^5 - x - 1/4
a, b = -9/8, 5/4
zs = find_zeros(f, (a,b))

p = plot(f, a, b; color=:black,
         linewidth=2, grid=false, legend=false, showaxis = false,
         background_color = :transparent)
plot!(p, zero, color=julia_colors[1])
scatter!(p, zs, zero.(zs), color=julia_colors[2:4], markersize=10)
save("src/assets/logo.png", p)
