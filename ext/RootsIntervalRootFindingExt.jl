module RootsIntervalRootFindingExt

using IntervalRootFinding
using Roots


function Roots.find_zeros(f, x0::IntervalRootFinding.Interval{T′}, M=IntervalRootFinding.Newton; kwargs...) where {T′}

    rts = IntervalRootFinding.roots(f, x0, M)
    T = float(T′)
    unique_roots = T[find_zero(f, (interval(r).lo, interval(r).hi)) for r ∈ rts if r.status == :unique]
    unknown = Interval{T}[interval(r) for r ∈ rts if r.status != :unique]

    (zeros = unique_roots, unknown=unknown)
end

Roots.find_zeros(f, x0::IntervalRootFinding.Interval, M::Roots.Newton; kwargs...) =
    Roots.find_zeros(f, x0, IntervalRootFinding.Newton)
Roots.find_zeros(f, x0::IntervalRootFinding.Interval, M::Roots.Bisection; kwargs...) =
    Roots.find_zeros(f, x0, IntervalRootFinding.Bisection)

end
