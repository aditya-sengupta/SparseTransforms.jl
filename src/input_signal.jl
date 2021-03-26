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
            wht[l+1] = s
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
        return new(n, k, locs, strengths, noise_sd, Dict{Int64,Float64}())
    end
end

Base.size(s::LazySignal) = length(s.signal_t)

function Base.getindex(s::LazySignal, i::Int)
    if haskey(s.signal_t, i)
        return s.signal_t[i]
    else
        m = dec_to_bin(i-1, s.n)
        sgns = (-1) .^ [m ⋅ dec_to_bin(l, s.n) for l in s.locs]
        v = s.strengths ⋅ sgns + rand(Normal(0.0, s.noise_sd))
        s.signal_t[i] = v
        return v
    end
end

function Base.getindex(s::Signal, r::AbstractUnitRange)
    return map(i -> s[i], r)
end

function Base.getindex(s::Signal, a::Array{<: Int})
    return map(i -> s[i], a)
end

function Base.getindex(s::Signal, a::Array{Array{<: Int}})
    return map(i -> s[i], a)
end

function Base.setindex!(s::LazySignal, v, i::Int)
    @warn("LazySignal is set based on transformed modes: explicitly setting a value is not supported.")
end

function Base.show(io::IO, s::LazySignal)
    println("Lazy signal with locations at ", s.locs)
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
Base.getindex(s::InputSignal, i::Int) = s.signal_t[i]
Base.setindex!(s::InputSignal, v, i::Int) = (s.signal_t[i] = v)

function get_random_sparse_signal(n::Int, K::Int, σ::Float64, minpower::Float64, maxpower::Float64; lazy::Bool=true)
    locs = sample(0:(2^n-1), K, replace=false)
    strengths = Float64.(rand(Uniform(minpower, maxpower), K)) .* (-1) .^ rand(Bool, K)
    if lazy
        return LazySignal(n, locs, strengths, σ)
    else
        return TestSignal(n, locs, strengths, σ)
    end
end

function get_random_delta_sparse_signal(n::Int, σ::Float64, delta::Float64, c::Int, minpower::Float64, maxpower::Float64)
    k = c * n ^ delta |> ceil |> Int64
    return get_random_sparse_signal(n, k, σ, minpower, maxpower)
end
