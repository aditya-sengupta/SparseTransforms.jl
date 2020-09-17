using Random
using Distributions

include("utils.jl")

struct InputSignal 
    n::Int64
    loc::Array{Int64,1}
    strengths::Array{Float64,1}
    noise_sd::Float64
    signal_w::Array{Float64,1}
    signal_t::Array{Float64,1}
    function InputSignal(n::Int64, loc::Array{Int64,1}, strengths::Array{Float64,1} = ones(length(loc)), noise_sd::Float64 = 0.0)
        N = 2^n
        noise = rand(Normal(0.0, noise_sd), N)
        wht = zeros(N)
        for (l, s) in zip(loc, strengths)
            wht[l] = s
        end
        signal_t = fwht(wht) + noise
        signal_w = fwht(signal_t) / N
        new(n, loc, strengths, noise_sd, signal_w, signal_t)
    end
end
