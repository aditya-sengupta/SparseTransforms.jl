using LinearAlgebra
using Pipe

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

function bin_to_dec(x::Array{Bool,1})
    """
    Takes in a binary array and returns the corresponding integer. 
    """
    return 2 .^(length(x)-1 : -1 : 0) ⋅ x
end

function dec_to_bin(x::Int64, num_bits::Int64)
    """
    Takes in an integer and returns the corresponding binary array.
    """
    @assert x < 2^num_bits "number of bits too small"
    u = split(lpad(string(x, base=2), num_bits, "0"), "")
    return map(y-> parse(Bool, y), u)
end

function binary_ints(m::Int64)
    """
    Returns a matrix where row 'i' is dec_to_bin(i, m), for i from 0 to 2 ** m - 1.
    """
    return hcat(map(x->dec_to_bin(x, m), 0:2^m-1)...)'
end

function sign_spright(x::Int64)
    """
    Sign function that matches the convention (footnote 2 on page 11):
    returns 0 for positive reals and 1 for negative reals.
    """
    return Int(x < 0)
end

function flip(x::Array{Bool,1})
    """
    Flip every bit in a boolean array.
    """
    return map(y -> !y, x)
end