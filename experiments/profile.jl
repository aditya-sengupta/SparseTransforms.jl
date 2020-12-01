using SparseTransforms
using ProgressMeter
using Plots

nrange = 4:20
num_runs = 5
methods = [:simple, :identity_like, :mle]
σ = 1e-2

function method_name(methods::Array{Symbol,1})
    if methods[3] == :nso
        return "non-sample-optimal"
    elseif methods[3] == :mle
        return "maximum likelihood"
    end
end

function bvals(n::Int64, btype::Symbol)
    if btype == :constant
        return 2
    elseif btype == :scaling
        if n < 8
            return 2
        else
            return n ÷ 4
        end
    end
end

function profile(methods::Array{Symbol,1}, btype::Symbol)
    brange = map(x -> bvals(x, btype), nrange)
    times = Float64[]
    samples = Float64[]
    println(brange)
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
    plot!(nrange, samples; label="experimental values", xlabel="log2(signal length)", ylabel="sample fraction", title="Sample ratio with $(method_name(methods)) decoding and $(string(btype)) sparsity")
end