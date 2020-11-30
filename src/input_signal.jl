using Random
using Distributions

include("utils.jl")

abstract type Signal end

struct TestSignal <: Signal
    n::Int64
    b::Int64
    loc::Array{Int64,1}
    strengths::Array{Float64,1}
    noise_sd::Float64
    signal_t::Array{Float64,1}

    function TestSignal(n::Int64, loc::Array{Int64,1}, strengths::Array{Float64,1} = ones(length(loc)), noise_sd::Float64 = 0.0)
        N = 2^n
        b = loc |> length |> log2 |> ceil |> Int64
        noise = rand(Normal(0.0, noise_sd), N)
        wht = zeros(N)
        for (l, s) in zip(loc, strengths)
            wht[l+1] = s # 1-indexing doesn't work here and I'm too lazy to use an offsetarray
        end
        signal_t = fwht(wht) + noise
        new(n, b, loc, strengths, noise_sd, signal_t)
    end
end

struct LazySignal <: Signal
    n::Int64
    b::Int64
    loc::Array{Int64,1}
    strengths::Array{Float64,1}
    noise_sd::Float64
    signal_t::Dict

    function LazySignal(n::Int64, loc::Array{Int64,1}, strengths::Array{Float64,1} = ones(length(loc)), noise_sd::Float64 = 0.0)
        N = 2^n
        b = loc |> length |> log2 |> ceil |> Int64
        new(n, b, loc, strengths, noise_sd, Dict())
    end
end

function get_subsignal(signal::Signal, ind::Int64)
    if isa(signal.signal_t, Array) || haskey(signal.signal_t, ind)
        return signal.signal_t[ind]
    else
        m = dec_to_bin(ind, signal.n)
        sgns = (-1) .^ [m ⋅ dec_to_bin(l, signal.n) for l in signal.loc]
        val = (1 / sqrt(2 ^ signal.n)) * (signal.strengths ⋅ sgns) + rand(Normal(0.0, signal.noise_sd))
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