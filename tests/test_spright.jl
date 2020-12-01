using SparseTransforms
using Test

function test_spright()
    println("Testing input signal constructor...")
    n, b = 4, 2
    signal = TestSignal(4, [4, 6, 10, 15])
    signal = TestSignal(4, [4, 6, 10, 15], Float64.([2, 4, 1, 1]))
    @test signal.n == n

    println("Testing query methods...")
    b = get_b(signal)
    Ms = get_Ms(n, b)
    D = get_D(n)
    whts_and_inds = map(M -> compute_delayed_subtransform(signal, M, D, fwht), Ms)

    println("Testing reconstruction...")
    bin_cardinalities = map(M -> bin_cardinality(signal, M, D), Ms)

    println("Testing full transform...")
    methods = [:simple, :nso, :nso, :none]
    spright_wht, utilization = spright(signal, methods; verbose=true, report=true)
    println("Used $(utilization / 2^n) of all time samples")

    println("Checking transform result...")
    true_wht = Dict(l => s for (l, s) in zip(signal.locs, signal.strengths))
    @test length(signal.locs) == length(spright_wht)
    for (i, x) in spright_wht
        @test true_wht[i] ≈ x
    end
    println("Transform is correct!")
end
