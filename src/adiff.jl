using PowerSeries
export D, D2

function D(f::Function, k::Int)
    @assert k > 0

    function(x)
        x = series(tuple(x, 1.0, ntuple(k-1, x->0.0)...)...)
        factorial(k) * f(x).(symbol(pop!(names(x)))) # or use polydir!
    end
end

D(f::Function) = D(f, 1)
D2(f::Function) = D(f, 2)
