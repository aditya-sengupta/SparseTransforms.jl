using Random
using Distributions

include("utils.jl")

abstract type Signal end

struct TestSignal <: Signal
    n::Int64
    loc::Array{Int64,1}
    strengths::Array{Float64,1}
    noise_sd::Float64
    signal_w::Array{Float64,1}
    signal_t::Array{Float64,1}

    function TestSignal(n::Int64, loc::Array{Int64,1}, strengths::Array{Float64,1} = ones(length(loc)), noise_sd::Float64 = 0.0)
        N = 2^n
        noise = rand(Normal(0.0, noise_sd), N)
        wht = zeros(N)
        for (l, s) in zip(loc, strengths)
            wht[l + 1] = s # 1-indexing doesn't work here and I'm too lazy to use an offsetarray
        end
        signal_t = fwht(wht) + noise
        signal_w = fwht(signal_t) / N
        new(n, loc, strengths, noise_sd, signal_w, signal_t)
    end
end

struct InputSignal <: Signal
    n::Int64
    signal_t::Array{Float64,1}

    function InputSignal(signal_t::Array{Float64,1})
        n = length(signal_t) |> log2 |> ceil |> Int64
        new(n, signal_t)
    end
end