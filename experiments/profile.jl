using SparseTransforms
using ProgressMeter

nrange = 4:2:20
num_runs = 3
methods = [:simple, :nso, :nso]
times = Float64[]
samples = Float64[]
σ = 1e-2

@showprogress for n = nrange
    t = 0
    s = 0
    b = 2
    println(n)
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
