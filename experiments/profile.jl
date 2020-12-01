using SparseTransforms
using ProgressMeter
using Plots

nrange = 4:20
num_runs = 5
methods = [:simple, :nso, :nso]
times = Float64[]
samples = Float64[]
σ = 1e-2

function bvals(n::Int64, btype::Symbol)
    if btype == :const
        return 2
    elseif btype == :scale
        if n < 8
            return n ÷ 2
        else
            return n ÷ 4
        end
    end
end

function profile(methods::Array{Symbol,1}, btype::Symbol)
    brange = map(x -> bvals(x, btype), nrange)
    @showprogress for (n, b) in zip(nrange, brange)
        t = 0
        s = 0
        for i = 1:num_runs
            locs = sample(0:2^n-1, 2^b, replace=false)
            strengths = Float64.(rand(Uniform(0.1, 10), 2^b)) .* (-1) .^ rand(Bool, 2^b)
            signal = LazySignal(n, locs, strengths, σ)
            res = @timed spright(signal, methods; report=true)
            s += res.value[2]
            t += res.time
        end
        push!(samples, s / (num_runs * 2 ^ n))
        push!(times, t / num_runs)
    end
    return samples, times
end

function plot_samples(methods::Array{Symbol,1}, btype::Symbol)
    brange = map(x -> bvals(x, btype), nrange)
    samples, times = profile(methods, btype)
    plot()
    plot!(nrange, (2 .^ brange) .* nrange ./ (2 .^ nrange); label="asymptotic trend")
    plot!(nrange, samples; label="experimental values", xlabel="log2(signal length)", ylabel="sample fraction", title="Sample efficiency with $(string.(methods)) and sparsity $(string(btype))")
end