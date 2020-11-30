using SparseTransforms
using Test

function test_signal()
    println("Testing signal retrieval")
    n, b = 4, 2
    signal = TestSignal(4, [4, 6, 10, 15], Float64.([2, 4, 1, 1]))

    println("Testing index retrieval")
    for i = 1:2^n
        @test get_subsignal(signal, i) == signal.signal_t[i]
    end

    subsample_check = collect(1:n:2^n) # or whatever subsample you want
    @test get_subsignal(signal, subsample_check) == map(i -> signal.signal_t[i], subsample_check)

    lazy_signal = LazySignal(4, [4, 6, 10, 15], Float64.([2, 4, 1, 1]))

    println("Testing lazy signal")
    for i = 1:2^n
        @test get_subsignal(signal, i) == get_subsignal(lazy_signal, i)
    end
end