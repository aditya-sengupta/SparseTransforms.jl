using SparseTransforms
using StatsBase, Random, Distributions
using Test

Random.seed!(1234)

function test_random()
    n = 8
    b = 4
    locs = sample(1:2^n, 2^b, replace=false)
    strengths = rand(Uniform(-10, 10), 2^b)
    signal = TestSignal(n, locs, strengths)

    methods = [:simple, :nso, :nso]
    spright_wht = spright(signal, methods; verbose=true)

    println("SPRIGHT result: ", spright_wht)
    for (l, s) in zip(signal.loc, signal.strengths)
        @test spright_wht[l] ≈ s
    end
end