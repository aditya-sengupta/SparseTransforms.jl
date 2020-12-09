using SparseTransforms
using Test

function test_signal()
    println("Testing signal retrieval")
    n, k = 8, 3
    locs = sample(0:2^n-1, 2^k, replace=false)
    strengths = Float64.(rand(Uniform(0.1, 10), 2^k)) .* (-1) .^ rand(Bool, 2^k)
    signal = TestSignal(n, locs, strengths)

    println("Testing index retrieval")
    for i = 1:2^n
        @test get_subsignal(signal, i) == signal.signal_t[i]
    end

    subsample_check = collect(1:n:2^n) # or whatever subsample you want
    @test get_subsignal(signal, subsample_check) == map(i -> signal.signal_t[i], subsample_check)

    lazy_signal = LazySignal(n, locs, strengths)

    println("Testing lazy signal")
    for i = 1:2^n
        @test get_subsignal(signal, i) ≈ get_subsignal(lazy_signal, i)
    end

    println("Testing multiple index subsampling")
    for _ = 1:20
        idxs = sample(1:2^n-1, 2^k, replace=false)
        @test get_subsignal(signal, idxs) ≈ get_subsignal(lazy_signal, idxs)
        # println(get_subsignal(signal, idxs))
        # println(get_subsignal(lazy_signal, idxs))
        # println("---------")
    end
end
