using SPRIGHT
using Test

println("Testing input signal constructor")
n, b = 4, 2
signal = InputSignal(4, [4, 6, 10, 15])
signal = InputSignal(4, [4, 6, 10, 15], Float64.([2, 4, 1, 1]))
@test signal.n == n
for (l, s) in zip(signal.loc, signal.strengths)
    @test signal.signal_w[l + 1] == s
end

println("Testing query methods...")
b = get_b(signal)
Ms = get_Ms(n, b)
D = get_D(n)
whts_and_inds = map(M -> compute_delayed_wht(signal, M, D), Ms)

println("Testing reconstruction...")
bin_cardinalities = map(M -> bin_cardinality(signal, M, D), Ms)

println("Testing full transform...")
methods = [:simple, :identity_like, :noiseless]
spright_wht = transform(signal, methods; verbose=true)

println("Checking transform result...")
@assert all(spright_wht .≈ signal.signal_w)
println("Transform is correct!")
