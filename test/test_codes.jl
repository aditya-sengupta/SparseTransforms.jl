using SparseTransforms
using Test

function test_codes()
    println("Testing LDPC code generation correctness")
    @test do_test(20, 25, 100)
end
