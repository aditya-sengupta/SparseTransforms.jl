using SparseTransforms
using Test

function test_query()
    n = 4
    num_delays = 6
    locs = [4, 6, 10, 15]
    signal = TestSignal(n, locs)

    # this also implicitly does type checks

    println("Testing all the query methods")
    for method_name in all_methods["query"]
        local b = get_b(signal; method=method_name)
        local Ms = get_Ms(n, b; method=method_name)
        @test b < n
        for M in Ms
            @test size(M) == (n, b)
        end
    end

    println("Testing the delays methods")
    for method_name in all_methods["delays"]
        local D = get_D(n; method=method_name, num_delays=num_delays)
        # @ test size(D, 1) > n TODO
        @test size(D, 2) == n
    end
end
