## add in simple derivative operator using PowerSeries package
## D(f), D(f,k) returns the derivative to power k (k=1 is default)
function D(f::Function, k::Int=1)
    @assert k > 0

    function(x)
        y = series(tuple(x, 1.0, ntuple(k-1, z->0.0)...)...)
        factorial(k) * f(y).(symbol(pop!(names(y)))) # or use polyder!
    end
end

D2(f::Function) = D(f, 2)
