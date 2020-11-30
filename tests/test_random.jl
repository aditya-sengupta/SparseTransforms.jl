using SparseTransforms
using StatsBase, Random, Distributions
using Test

Random.seed!(1234)

function test_random()
    n = 4
    b = 2
    σ = 1e-2
    locs = sample(0:2^n-1, 2^b, replace=false)
    strengths = Float64.(rand(Uniform(0.1, 10), 2^b)) .* (-1) .^ rand(Bool, 2^b)
    signal = TestSignal(n, locs, strengths, σ)

    println("True locations: ", locs)
    println("True strengths: ", strengths)

    other_b = get_b(signal; method=:simple)
    Ms = get_Ms(n, b; method=:simple)
    for loc in locs
        println(loc)
        for (i, M) in enumerate(Ms)
            exp = expected_bin(dec_to_bin(loc, n), M)
            println("\t$i: $exp")
        end
        println()
    end

    methods = [:simple, :nso, :nso]
    spright_wht, used_size = spright(signal, methods; verbose=true, report=true)

    println("SPRIGHT result: ", spright_wht)
    println("Used $(used_size / 2^n) of all time samples")
    @test length(signal.loc) == length(spright_wht)
    for (l, s) in zip(signal.loc, signal.strengths)
        @test isapprox(spright_wht[l], s, atol=5*σ)
    end
end
