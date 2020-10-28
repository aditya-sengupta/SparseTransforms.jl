using SPRIGHT
using Test

println("Testing input signal constructor")
n, b = 4, 2
signal = InputSignal(4, [4, 6, 10, 15])
signal = InputSignal(4, [4, 6, 10, 15], Float64.([2, 4, 1, 1]))
@test signal.n == n
for (l, s) in zip(signal.loc, signal.strengths)
    @test signal.signal_w[l] == s
end

println("Testing query methods")
b = get_b(signal)
Ms = get_Ms(n, b)
D = get_D(n)


