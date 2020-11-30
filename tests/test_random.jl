using SparseTransforms
using StatsBase, Random, Distributions
using Test

Random.seed!(1234)

function test_random()
    n = 16
    k = 2
    σ = 1e-2
    locs = sample(0:2^n-1, 2^k, replace=false)
    strengths = Float64.(rand(Uniform(5, 10), 2^k)) .* (-1) .^ rand(Bool, 2^k)
    signal = LazySignal(n, locs, strengths, σ)

    println("True locations: ", locs)
    println("True strengths: ", strengths)

    b = get_b(signal; method=:simple)
    Ms = get_Ms(n, b; method=:simple)
    contents = Dict()
    for loc in locs
        println(loc)
        for (i, M) in enumerate(Ms)
            j = expected_bin(dec_to_bin(loc, n), M)
            if haskey(contents, (i,j))
                append!(contents[(i,j)], loc)
            else
                contents[(i,j)] = [loc]
            end
            println("\t$i: $j")
        end
        println()
    end
    println("true contents:")
    for (k,v) in contents
        println("\t$k: $v")
    end

    methods = [:simple, :nso, :nso]
    spright_wht, used_size = spright(signal, methods; verbose=true, report=true)

    println("SPRIGHT result: ", spright_wht)
    println("Used $(used_size / 2^n) of all time samples")
    @test length(signal.loc) == length(spright_wht)
    for (l, s) in zip(signal.loc, signal.strengths)
        println(spright_wht[l] / s)
        # @test isapprox(spright_wht[l], s, atol=5*σ)
    end
end
