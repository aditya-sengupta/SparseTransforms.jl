using Random
using Distributions

include("utils.jl")

abstract type Signal end

struct TestSignal <: Signal
    n::Int64
    k::Int64
    locs::Array{Int64,1}
    strengths::Array{Float64,1}
    noise_sd::Float64
    signal_t::Array{Float64,1}

    function TestSignal(n::Int64, locs::Array{Int64,1}, strengths::Array{Float64,1} = ones(length(locs)), noise_sd::Float64 = 0.0)
        N = 2^n
        k = locs |> length |> log2 |> ceil |> Int64
        noise = rand(Normal(0.0, noise_sd), N)
        wht = zeros(N)
        for (l, s) in zip(locs, strengths)
            wht[l+1] = s # 1-indexing doesn't work here and I'm too lazy to use an offsetarray
        end
        signal_t = fwht(wht) + noise
        new(n, k, locs, strengths, noise_sd, signal_t)
    end
end

struct LazySignal <: Signal
    n::Int64
    k::Int64
    locs::Array{Int64,1}
    strengths::Array{Float64,1}
    noise_sd::Float64
    signal_t::Dict

    function LazySignal(n::Int64, locs::Array{Int64,1}, strengths::Array{Float64,1} = ones(length(locs)), noise_sd::Float64 = 0.0)
        N = 2^n
        k = locs |> length |> log2 |> ceil |> Int64
        new(n, k, locs, strengths, noise_sd, Dict())
    end
end

function get_subsignal(signal::Signal, ind::Int64)
    if isa(signal.signal_t, Array) || haskey(signal.signal_t, ind)
        return signal.signal_t[ind]
    else
        m = dec_to_bin(ind-1, signal.n)
        sgns = (-1) .^ [m ⋅ dec_to_bin(l, signal.n) for l in signal.locs]
        # val = (1 / sqrt(2 ^ signal.n)) * (signal.strengths ⋅ sgns) + rand(Normal(0.0, signal.noise_sd))
        val = signal.strengths ⋅ sgns + rand(Normal(0.0, signal.noise_sd))
        signal.signal_t[ind] = val
        return val
    end
end

function get_subsignal(signal::Signal, inds::Array{Int64,1})
    return map(x -> get_subsignal(signal, x), inds)
end

struct InputSignal <: Signal
    n::Int64
    signal_t::Array{Float64,1}
    noise_sd::Float64

    function InputSignal(signal_t::Array{Float64,1}, noise_sd::Float64)
        n = length(signal_t) |> log2 |> ceil |> Int64
        new(n, signal_t, noise_sd)
    end
end

function get_random_sparse_signal(n::Int64, k::Int64, σ::Float64; lazy::Bool=true)
    locs = sample(0:2^n-1, 2^k, replace=false)
    strengths = Float64.(rand(Uniform(5, 10), 2^k)) .* (-1) .^ rand(Bool, 2^k)
    if lazy
        return LazySignal(n, locs, strengths, σ)
    else
        return TestSignal(n, locs, strengths, σ)
    end
end
