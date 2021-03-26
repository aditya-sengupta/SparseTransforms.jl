"""
Methods to profile the SPRIGHT algorithm and this implementation of it.
"""

using SparseTransforms
using ProgressMeter
using Plots
using Profile

nrange = 4:20
num_runs = 5
methods_list = [:simple, :nso, :nso, :none]
σ = 1e-2
deltas = 0.05:0.05:0.3 |> collect

function method_name(methods_list::Array{Symbol,1})
    if methods_list[3] == :nso
        return "non-sample-optimal"
    elseif methods_list[3] == :mle
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

function profile_fixed_b(methods_list::Array{Symbol,1}, btype::Symbol)
    brange = map(x -> bvals(x, btype), nrange)
    times = Float64[]
    samples = Float64[]
    @showprogress for (n, b) in zip(nrange, brange)
        t = 0
        s = 0
        for i = 1:num_runs
            # locs = sample(0:2^n-1, 2^b, replace=false)
            # strengths = Float64.(rand(Uniform(0.1, 10), 2^b)) .* (-1) .^ rand(Bool, 2^b)
            # signal = LazySignal(n, locs, strengths, σ)
            signal = get_random_sparse_signal(n, 2^b, σ, 0.1, 10)
            res = @timed spright(signal, methods_list; report=true)
            s += res.value[2]
            t += res.time
        end
        push!(samples, s / (num_runs * 2 ^ n))
        push!(times, t / num_runs)
    end
    return samples, times
end

function profile_delta_sparse(methods_list::Array{Symbol,1}, deltas::Array{Float64,1})
    times_arr = Any[]
    samples_arr = Any[]
    @showprogress for δ in deltas
        times = Float64[]
        samples = Float64[]
        for n in nrange
            s, t = 0, 0
            for i = 1:num_runs
                signal = get_random_delta_sparse_signal(n, σ, δ, 1, 0.5, 10.0)
                res = @timed spright(signal, methods_list; report=true)
                s += res.value[2]
                t += res.time
            end
            push!(samples, s / (num_runs * 2 ^ n))
            push!(times, t / num_runs)
        end
        push!(times_arr, times)
        push!(samples_arr, samples)
    end
    return samples_arr, times_arr
end

function plot_samples(methods_list::Array{Symbol,1}, btype::Symbol)
    brange = map(x -> bvals(x, btype), nrange)
    samples, times = profile(methods_list, btype)
    plot()
    plot!(nrange, (2 .^ brange) .* nrange ./ (2 .^ nrange); label="asymptotic trend")
    plot!(nrange, samples; label="experimental values", xlabel="log2(signal length)", ylabel="sample fraction", title="Sample ratio with $(method_name(methods_list)) decoding and $(string(btype)) sparsity")
end

function plot_delta_profile()
    samples, times = profile_delta_sparse(methods_list, deltas)
    plot()
    plot!(nrange, samples; label=hcat(deltas...), xlabel="log2(signal length)", ylabel="sample fraction", title="Sample ratio with sparsity O(N^delta)")
    brange = map(x -> 2^(-0.7x) * x^2 / 12.5, nrange)
    plot!(nrange, brange; label="asymptotic trend")
end

function plot_delta_times()
    samples, times = profile_delta_sparse(methods_list, deltas)
    plot()
    plot!(nrange, times; label=hcat(deltas...), xlabel="log2(signal length)", ylabel="runtime", title="Runtime with sparsity O(N^delta)")
    brange = map(x -> 2^(0.3x) * x^3, nrange)
    # println(times[1,1])
    # println(times[1][1])
    plot!(nrange, brange * times[1][1] / brange[1] / 5; label="asymptotic trend")
end

function plot_delta_times_m_diff()
    samples, times = profile_delta_sparse(methods_list, [0.3])
    plot()
    plot!(nrange, times; label="identity subsampling", xlabel="log2(signal length)", ylabel="runtime", title="Runtime with sparsity O(N^delta)")
    methods_list[1] = :random
    samples2, times2 = profile_delta_sparse(methods_list, [0.3])
    plot!(nrange, times2; label="random subsampling")
    brange = map(x -> 2^(0.3x) * x^3, nrange)
    plot!(nrange, brange * times[1][1] / brange[1] / 5; label="asymptotic trend")
end
