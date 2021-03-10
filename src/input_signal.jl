using Random
using Distributions

include("utils.jl")

abstract type Signal <: AbstractArray{Float64,1} end # maybe extend to N-d SPRIGHT? 
Base.IndexStyle(::Signal) = IndexLinear()

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

Base.size(s::TestSignal) = size(s.signal_t)
Base.getindex(s::TestSignal, i::Int64) = s.signal_t[i]
Base.setindex!(s::TestSignal, v, i::Int64) = (s.signal_t[i] = v)

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

Base.size(s::LazySignal) = length(s.signal_t)
Base.getindex(s::LazySignal, i::Int64) = s.signal_t[i]
Base.setindex!(s::LazySignal, v, i::Int64) = (s.signal_t[i] = v)

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

Base.size(s::InputSignal) = size(s.signal_t)
Base.getindex(s::InputSignal, i::Int64) = s.signal_t[i]
Base.setindex!(s::InputSignal, v, i::Int64) = (s.signal_t[i] = v)

function get_random_sparse_signal(n::Int64, K::Int64, σ::Float64, minpower::Float64, maxpower::Float64; lazy::Bool=true)
    locs = sample(0:2^n-1, K, replace=false)
    strengths = Float64.(rand(Uniform(minpower, maxpower), K)) .* (-1) .^ rand(Bool, K)
    if lazy
        return LazySignal(n, locs, strengths, σ)
    else
        return TestSignal(n, locs, strengths, σ)
    end
end

function get_random_delta_sparse_signal(n::Int64, σ::Float64, delta::Float64, c::Int64, minpower::Float64, maxpower::Float64)
    k = c * n ^ delta |> ceil |> Int64
    return get_random_sparse_signal(n, k, σ, minpower, maxpower)
end
