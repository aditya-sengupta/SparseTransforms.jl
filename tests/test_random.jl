using SparseTransforms
using StatsBase, Random, Distributions
using Test

<<<<<<< HEAD
Random.seed!(1234)

function test_random()
    n = 8
    b = 4
    locs = sample(1:2^n, 2^b, replace=false)
    strengths = Float64.(rand(1:10, 2^b)) .* (-1) .^ rand(Bool, 2^b)
    signal = TestSignal(n, locs, strengths)
=======
function test_random()
    Random.seed!(1234)
    n = 16
    b = 2
    σ = 1e-2
    locs = sample(1:2^n, 2^b, replace=false)
    strengths = Float64.(rand(Uniform(0.1, 10), 2^b)) .* (-1) .^ rand(Bool, 2^b)
    signal = TestSignal(n, locs, strengths, σ)
>>>>>>> 902cf986dae27c3786e7058cd42477592fae94de

    println("True locations: ", locs)
    println("True strengths: ", strengths)

<<<<<<< HEAD
    methods = [:simple, :nso, :nso]
    spright_wht = spright(signal, methods; verbose=true)

    println("True result: ", sort(Dict(l => s for (l, s) in zip(locs, strengths))))
    println("SPRIGHT result: ", sort(spright_wht))
    for (l, s) in zip(signal.loc, signal.strengths)
        @test spright_wht[l] ≈ s
    end
end

function get_to_fail()
    while true
        test_random()
=======
    other_b = get_b(signal; method=:simple)
    Ms = get_Ms(n, other_b; method=:simple)
    for loc in locs
        println(loc)
        for (i, M) in enumerate(Ms)
            exp = expected_bin(dec_to_bin(loc, n), M)
            println("\t$i: $exp")
        end
        println()
    end

    methods = [:simple, :nso, :nso]
    spright_wht = spright(signal, methods; verbose=true)

    println("SPRIGHT result: ", spright_wht)
    @test length(signal.loc) == length(spright_wht)
    for (l, s) in zip(signal.loc, signal.strengths)
        @test isapprox(spright_wht[l], s, atol=5*σ)
>>>>>>> 902cf986dae27c3786e7058cd42477592fae94de
    end
end
