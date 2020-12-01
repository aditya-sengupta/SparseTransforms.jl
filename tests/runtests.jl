using SparseTransforms
using Test
include("test_query.jl")
include("test_reconstruct.jl")
include("test_spright.jl")

function tests()
    println("Testing SPRIGHT")
    test_query()
    test_reconstruct()
    test_spright()
end
