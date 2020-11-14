using SparseTransforms
using Test

println("Testing SPRIGHT")
include("test_query.jl")
include("test_reconstruct.jl")
include("test_spright.jl")
