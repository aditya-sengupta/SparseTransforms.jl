using SparseTransforms
using Test

function test_spright()
    println("Testing input signal constructor...")
    n, b = 4, 2
    signal = TestSignal(4, [4, 6, 10, 15])
    signal = TestSignal(4, [4, 6, 10, 15], Float64.([2, 4, 1, 1]))
    @test signal.n == n
    wht = fwht(signal.signal_t)
    for (l, s) in zip(signal.loc, signal.strengths)
        @test wht[l + 1] == s * 2^(signal.n)
    end

    println("Testing query methods...")
    b = get_b(signal)
    Ms = get_Ms(n, b)
    D = get_D(n)
    whts_and_inds = map(M -> compute_delayed_subtransform(signal, M, D, fwht), Ms)

    println("Testing reconstruction...")
    bin_cardinalities = map(M -> bin_cardinality(signal, M, D), Ms)

    println("Testing full transform...")
    methods = [:simple, :nso, :nso]
    spright_wht = spright(signal, methods; verbose=true)

    println("Checking transform result...")
    @test length(signal.loc) == length(spright_wht)
    for (l, s) in zip(signal.loc, signal.strengths)
        @test spright_wht[l] ≈ s
    end
    println("Transform is correct!")
end
