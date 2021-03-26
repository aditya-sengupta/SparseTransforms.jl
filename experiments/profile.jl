using SparseTransforms
using ProgressMeter
using Plots

nrange = 4:16
num_runs = 40
def_methods = [:simple, :nso, :nso, :none]
σ = 1.0
deltas = 0.05:0.05:0.3 |> collect

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

function profile_fixed_b(methods::Array{Symbol,1}, btype::Symbol)
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
            res = @timed spright(signal, methods; report=true)
            s += res.value[2]
            t += res.time
        end
        push!(samples, s / (num_runs * 2 ^ n))
        push!(times, t / num_runs)
    end
    return samples, times
end

function check_transform_correctness(target, estimated, σ)
    if length(target.locs) != length(estimated)
        return false
    end
    for (l, s) in zip(target.locs, target.strengths)
        if !haskey(estimated, l) || !isapprox(estimated[l], s, atol=2*σ)
            return false
        end
    end
    return true
end

function profile_delta_sparse(methods::Array{Symbol,1}, deltas::Array{Float64,1})
    times_arr = Any[]
    samples_arr = Any[]
    correct_arr = Any[]
    @showprogress for δ in deltas
        times = Float64[]
        samples = Float64[]
        correct = Float64[]
        for n in nrange
            s, t, c = 0, 0, 0
            for i = 1:num_runs
                signal = get_random_delta_sparse_signal(n, σ, δ, 1, 0.5, 10.0)
                res = @timed spright(signal, methods; report=true)
                s += res.value[2]
                t += res.time
                c += check_transform_correctness(signal, res.value[1], σ)
            end
            push!(samples, s / (num_runs * 2 ^ n))
            push!(times, t / num_runs)
            push!(correct, c / num_runs)
        end
        push!(times_arr, times)
        push!(samples_arr, samples)
        push!(correct_arr, correct)
    end
    return samples_arr, times_arr, correct_arr
end

function plot_samples(methods::Array{Symbol,1}, btype::Symbol)
    brange = map(x -> bvals(x, btype), nrange)
    samples, times = profile(methods, btype)
    plot()
    plot!(nrange, (2 .^ brange) .* nrange ./ (2 .^ nrange); label="asymptotic trend")
    plot!(nrange, samples; label="experimental values", xlabel="log2(signal length)", ylabel="sample fraction", title="Sample ratio with $(method_name(methods)) decoding and $(string(btype)) sparsity")
end

function plot_delta_profile()
    samples, times = profile_delta_sparse(def_methods, deltas)
    plot()
    plot!(nrange, samples; label=hcat(deltas...), xlabel="log2(signal length)", ylabel="sample fraction", title="Sample ratio with sparsity O(N^delta)")
    brange = map(x -> 2^(-0.7x) * x^2 / 12.5, nrange)
    plot!(nrange, brange; label="asymptotic trend")
end

function plot_delta_times()
    samples, times = profile_delta_sparse(def_methods, deltas)
    plot()
    plot!(nrange, times; label=hcat(deltas...), xlabel="log2(signal length)", ylabel="runtime", title="Runtime with sparsity O(N^delta)")
    brange = map(x -> 2^(0.3x) * x^3, nrange)
    # println(times[1,1])
    # println(times[1][1])
    plot!(nrange, brange * times[1][1] / brange[1] / 5; label="asymptotic trend")
end

function plot_delta_times_m_diff()
    samples, times = profile_delta_sparse(def_methods, [0.3])
    plot()
    plot!(nrange, samples; label="identity subsampling", xlabel="log2(signal length)", ylabel="Sample ratio", title="Sample ratio with sparsity O(N^0.3)")
    def_methods[1] = :random
    samples2, times2 = profile_delta_sparse(def_methods, [0.3])
    plot!(nrange, samples2; label="random subsampling")
    brange = map(x -> 2^(-0.7x) * x^2, nrange)
    plot!(nrange, brange / brange[1]; label="asymptotic trend")
end

function N_logb_N(b)
    return n -> 2^n * n^b
end

function K_logb_N(delta, b)
    return n -> 2^(delta * n) * n^b
end

function plot_method_comparison(plot_val)
    trial_ms = ([:simple, :identity_like, :noiseless, :none],
                [:simple, :identity_like, :mle, :none],
                [:simple, :nso, :nso, :none])
    delta = 0.3
    plot()
    samples_list = []
    times_list = []

    m = trial_ms[1]
    samples, times, corrects = profile_delta_sparse(m, [delta])
    append!(samples_list, samples)
    append!(times_list, times)
    plot_dict = Dict(:samples=>(samples, "Sample ratio"),
                     :times=>(times, "Runtime"),
                     :corrects=>(corrects, "Success rate"))
    y, ylab = plot_dict[plot_val]
    println(y)
    plot!(nrange, y; label=string(m[3]), xlabel="log2(signal length)", ylabel=ylab, title=ylab*" with sparsity O(N^$delta)")

    m = trial_ms[2]
    samples, times, corrects = profile_delta_sparse(m, [delta])
    append!(samples_list, samples)
    append!(times_list, times)
    plot_dict = Dict(:samples=>(samples, "Sample ratio"),
                     :times=>(times, "Runtime"),
                     :corrects=>(corrects, "Success rate"))
    y, ylab = plot_dict[plot_val]
    println(y)
    plot!(nrange, y; label=string(m[3]), xlabel="log2(signal length)", ylabel=ylab, title=ylab*" with sparsity O(N^$delta)")

    m = trial_ms[3]
    samples, times, corrects = profile_delta_sparse(m, [delta])
    append!(samples_list, samples)
    append!(times_list, times)
    plot_dict = Dict(:samples=>(samples, "Sample ratio"),
                     :times=>(times, "Runtime"),
                     :corrects=>(corrects, "Success rate"))
    y, ylab = plot_dict[plot_val]
    println(y)
    plot!(nrange, y; label=string(m[3]), xlabel="log2(signal length)", ylabel=ylab, title=ylab*" with sparsity O(N^$delta)")

    # if plot_val != :corrects
    #     param_lookup = Dict(:samples => (2, delta-1, samples_list[3][10]),
    #                         :times => (3, delta, times_list[3][10]))
    #     b, d, norm = param_lookup[plot_val]
    #     asymp = map(K_logb_N(d, b), nrange)
    #     plot!(nrange, asymp * norm / asymp[10]; label="asymptotic trend")
    # end
    #
    # gui()

    # for m in trial_ms
    #     println(m)
    #     s, t = profile_delta_sparse(m, [0.3])
    #     println(t)
    #     plot!(nrange, t; label=string(m[3]), xlabel="log2(signal length)", ylabel="runtime", title="Runtime with sparsity O(N^0.3)")
    #     # samples_trials.append(s)
    #     # times_trials.append(t)
    # end
end
