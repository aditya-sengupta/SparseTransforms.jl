using SparseTransforms
using Test
include("test_signal.jl")
include("test_query.jl")
include("test_reconstruct.jl")
include("test_spright.jl")
include("test_random.jl")
include("test_codes.jl")

function tests()
    println("Testing SPRIGHT")
    test_signal()
    test_query()
    test_reconstruct()
    test_spright()
    test_random()
    test_codes()
end

tests()