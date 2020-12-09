using SparseTransforms
using StatsBase, Random, Distributions
using Test

Random.seed!(1234)

function test_random()
    n = 8
    k = 2
    σ = 1e-5
    minpower, maxpower = 0.5, 10.0
    signal = get_random_sparse_signal(n, 2^k, σ, minpower, maxpower)

    println("True locations: ", signal.locs)
    println("True strengths: ", signal.strengths)

    b = get_b(signal; method=:simple)
    Ms = get_Ms(n, b; method=:simple)
    contents = Dict()
    for loc in signal.locs
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

    methods = [:simple, :nso, :nso, :none]
    # methods = [:simple, :random, :mle, :none]
    spright_wht, used_size = spright(signal, methods; verbose=true, report=true)

    println("SPRIGHT result: ", spright_wht)
    println("Used $(used_size / 2^n) of all time samples")
    @test length(signal.locs) == length(spright_wht)
    for (l, s) in zip(signal.locs, signal.strengths)
        # println(spright_wht[l] / s)
        @test isapprox(spright_wht[l], s, atol=5*σ)
    end
end
