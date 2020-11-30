using SparseTransforms
using ProgressMeter
using Distributed

nrange = [4, 8, 16]
num_runs = 10
methods = [:simple, :nso, :nso]
times = Float64[]
samples = Float64[]
σ = 1e-2

@showprogress @distributed for n = nrange
    t = 0
    s = 0
    b = n ÷ 2
    
    for i = 1:num_runs
        locs = sample(0:2^n-1, 2^b, replace=false)
        strengths = Float64.(rand(Uniform(0.1, 10), 2^b)) .* (-1) .^ rand(Bool, 2^b)
        signal = TestSignal(n, locs, strengths, σ)
        res = @timed spright(signal, methods; report=true)
        println(res)
        s += res.value[2]
        t += res.time
    end
    push!(samples, s / num_runs)
    push!(times, t / num_runs)
end
