"""
Methods for the query generator: specifically, to

1. generate sparsity coefficients b and subsampling matrices M
2. get the indices of a signal subsample
3. compute a subsampled and delayed Walsh-Hadamard transform.
"""

# potential TODO figure out how to load in parts of the signal bit by bit
# see how Amirali did it

include("input_signal.jl")
using StatsBase
using Random

"""
A semi-arbitrary fixed choice of the sparsity coefficient.
See get_b for full signature.
"""
function get_b_simple(signal::Signal)
    # I don't know why, but changing this to signal.n \div 3 makes it break
    if signal.n == 4
        return 2
    elseif (signal.n == 8) || (signal.n == 16)
        return 4
    elseif (signal.n < 16)
        return signal.n ÷ 2
    else
        return signal.n ÷ 4
    end
end

get_b_lookup = Dict(
    :simple => get_b_simple,
    :random => get_b_simple
)

"""
Get the sparsity coefficient for the signal.

Arguments
---------
signal : TestSignal
The signal whose WHT we want.

method : Symbol
The method to use. All methods referenced must use this signature (minus "method".)

Returns
-------
b : Int64
The sparsity coefficient.
"""
function get_b(signal::Signal; method=:simple)::Int64
    return get_b_lookup[method](signal)
end

function get_partial_I(n, b, k)
    M = falses(n, b)
    for j in 1:b
        M[k+j, j] = 1
    end
    return M
end

"""
A semi-arbitrary fixed choice of the subsampling matrices. See get_Ms for full signature.
"""
function get_Ms_simple(n::Int64, b::Int64)
    # @assert n % b == 0 "n must be exactly divisible by b"
    num_to_get = Int64(ceil(n / b))

    Ms = BitArray[]
    push!(Ms, get_partial_I(n, b, n-b))
    for i in num_to_get-2:-1:0
        push!(Ms, get_partial_I(n, b, b*i))
    end
    return Ms
end

function get_Ms_random(n::Int64, b::Int64)
    num_to_get = n ÷ b
    Ms = BitArray[]

    for i in 1:num_to_get
        push!(Ms, bitrand(n, b))
    end
    return Ms
end

get_Ms_lookup = Dict(
    :simple => get_Ms_simple,
    :random => get_Ms_random
)

function get_Ms(n::Int64, b::Int64; method=:simple)::Array{BitArray,1}
    return get_Ms_lookup[method](n, b)
end

"""
Gets the delays matrix [0; I], of dimension (n+1, n). See get_D for full signature.
"""
function get_D_identity_like(n::Int64; kwargs...)
    return BitArray(vcat(zeros(1, n), Matrix{Bool}(I, n, n)))
end

"""
Gets a random delays matrix of dimension (num_delays, n). See get_D for full signature.
"""
function get_D_random(n::Int64; kwargs...)
    num_delays = kwargs[:num_delays]
    choices = sample(0:2^n-1, num_delays, replace=false)
    return @pipe dec_to_bin.(choices, n) |> hcat(_...) |> transpose |> BitArray
end

"""
Get a repetition code based (NSO-SPRIGHT) delays matrix. See get_D for full signature.
Not sure of correctness based on the paper: this is not (num_delays, n).
"""
function get_D_nso(n::Int64; kwargs...)::BitArray{2}
    num_delays = kwargs[:num_delays]
    p1 = num_delays ÷ n
    random_offsets = get_D_random(n; num_delays=p1)
    D = BitArray(undef, 0, n)
    identity_like = get_D_identity_like(n)
    for row in eachrow(random_offsets)
        modulated_offsets = @pipe [(row .⊻ r) for r in eachrow(identity_like)] |> hcat(_...) |> transpose
        D = vcat(D, modulated_offsets)
    end
    return D
end

function get_generator_matrix(codelen::Int64, n::Int64, code::Symbol)
    return zeros((codelen, n))
end

"""
Get a channel code based (SO-SPRIGHT) delays matrix. See get_D for full signature.
Returns a matrix of dimension (p1+p2+p3, n).
"""
function get_D_so(n::Int64; kwargs...)::BitArray{2}
    p1, p2, p3 = kwargs[:num_delays]
    code = kwargs[:code]
    D1 = get_D_random(n, num_delays=p1)
    D2 = zeros((p2, n))
    D3 = get_generator_matrix(p3, n, code)
    return vcat(D1, D2, D3)
end

get_D_lookup = Dict(
    :identity_like => get_D_identity_like,
    :random => get_D_random,
    :nso => get_D_nso,
    :so => get_D_so
)

"""
Delay generator: gets a delays matrix.

Arguments
---------
n : Int64

number of bits: log2 of the signal length.

Returns
-------
D : num_delays×n Array{Bool,2}

The delays matrix; if num_delays is not specified in kwargs, see the relevant sub-function for a default.
"""
function get_D(n::Int64; method=:identity_like, kwargs...)
    return get_D_lookup[method](n; kwargs...)
end

"""
Query generator: creates indices for signal subsamples.
Arguments
---------
M : n×b Array{Bool,2}
The subsampling matrix.

d : n-element Array{Bool,1}
The subsampling offset.

Returns
-------
indices : B-element Array{Int64,1}
The (decimal) subsample indices.
"""
function subsample_indices(M, d)
    L = binary_ints(size(M, 2))
    inds_binary = Bool.((M*L .⊻ d) .% 2)
    return @pipe inds_binary |> eachcol |> collect |> bin_to_dec.(_) .+ 1
end

"""
Creates random delays, subsamples according to M and the random delays,
and returns the subsample WHT along with the delays.

Arguments
---------
signal : Signal
The signal to subsample, delay, and compute the WHT of.

M : n×b Array{Bool,2}
The subsampling matrix.

num_delays : Int64
The number of delays to apply; or, the number of rows in the delays matrix.

force_identity_like : Bool
Whether to make D = [0; I] like in the noiseless case; for debugging.
"""
function compute_delayed_wht(signal::Signal, M, D)
    return compute_delayed_subtransform(signal, M, D, fwht)
end

function compute_delayed_fft(signal::Signal, M, D)
    return compute_delayed_subtransform(signal, M, D, fft)
end

function compute_delayed_subtransform(signal::Signal, M, D, transform::Function)
    inds = map(d -> subsample_indices(M, d), D |> eachrow |> collect)
    used_inds = reduce(union, inds)
    samples_to_transform = map(x -> signal[x], inds)
    return hcat(transform.(samples_to_transform)...), used_inds
end
