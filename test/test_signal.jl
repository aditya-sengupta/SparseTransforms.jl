using SparseTransforms
using Test

function test_signal()
    println("Testing signal retrieval")
    n, k = 8, 3
    locs = sample(0:2^n-1, 2^k, replace=false)
    strengths = Float64.(rand(Uniform(0.1, 10), 2^k)) .* (-1) .^ rand(Bool, 2^k)
    signal = TestSignal(n, locs, strengths)

    lazy_signal = LazySignal(n, locs, strengths)

    println("Testing lazy signal")
    for i = 1:2^n
        @test signal[i] ≈ lazy_signal[i]
    end

    println("Testing multiple index subsampling")
    for _ = 1:20
        idxs = sample(1:2^n-1, 2^k, replace=false)
        @test signal[idxs] ≈ lazy_signal[idxs]
    end
end
