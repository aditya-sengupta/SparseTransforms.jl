"""
Linear channel coding methods for frequency detection; specifically to

1. create generator matrices
2. BSC decode the results
"""

get_code_lookup = Dict()

"""
Estimate the frequency from a singleton given an offset coding scheme and generator matrix G.

Arguments
---------
recvd_codeword : BitArray{1}
The sampled vector from which we are decoding

G : BitArray{2}
The generator matrix fo the code

code : Symbol
The code to use. All methods referenced must use this signature (minus "code".)

Returns
-------
k : BitArray{1}
The estimated frequency.
"""
function decode_with(recvd_codeword::BitArray{1}, G::BitArray{2}, code::Symbol)::BitArray{1}
    return get_code_lookup[code](recvd_codeword, G)
end
