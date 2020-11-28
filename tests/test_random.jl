using SparseTransforms
using StatsBase, Random, Distributions
using Test

Random.seed!(1234)

function test_random()
    n = 4
    b = 2
    locs = sample(1:2^n, 2^b, replace=false)
    strengths = Float64.(rand(1:10, 2^b)) .* (-1) .^ rand(Bool, 2^b)
    signal = TestSignal(n, locs, strengths)

    println("True locations: ", locs)
    println("True strengths: ", strengths)

    methods = [:simple, :nso, :nso]
    spright_wht = spright(signal, methods; verbose=true)

    println("SPRIGHT result: ", spright_wht)
    for (l, s) in zip(signal.loc, signal.strengths)
        @test spright_wht[l] ≈ s
    end
end