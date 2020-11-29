using LinearAlgebra
using FFTW
using Pipe.Pipe

"""
Cooley-Tukey recursive Fast Walsh-Hadamard Transform.
"""
function fwht(x::Array{Float64,1})
    N = length(x)
    if N == 1
        return x
    else
        X_even = fwht(x[1 : N÷2])
        X_odd = fwht(x[N÷2+1 : N])
        return cat(X_even + X_odd, X_even - X_odd, dims=1)
    end
end

"""
Takes in a binary array and returns the corresponding integer.
"""
function bin_to_dec(x::AbstractArray{Bool})
    return 2 .^(length(x)-1 : -1 : 0) ⋅ x
end

"""
Takes in an integer and returns the corresponding binary array.
"""
function dec_to_bin(x::Int64, num_bits::Int64)
    @assert x < 2^num_bits "number of bits too small"
    return @pipe string(x, base=2) |> lpad(_, num_bits, "0") |> split(_, "") |> parse.(Bool, _)
end

"""
Returns a matrix where row 'i' is dec_to_bin(i, m),
for i from 0 to 2 ** m - 1.
"""
function binary_ints(m::Int64)
    return @pipe 0:2^m-1 |> dec_to_bin.(_, m) |> hcat(_...)
end

"""
Sign function that matches the convention (footnote 2 on page 11):
returns 0 for positive reals and 1 for negative reals.
"""
function sign_spright(x)
    return Int64(x < 0)
end

function expected_bin(k::BitArray{1}, M::BitArray{2})
    return @pipe M'k |> Bool.(_) |> bin_to_dec(_) + 1
end
