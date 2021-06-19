using Random

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

#LDPC Code Based on https://github.com/christianpeel/LDPCsim.jl/blob/master/src/LDPCsim.jl

function generate_LDPC_adjacency_matrix(M::Int64, N::Int64)
    # Create LDPC adjacency matrix of the Tanner graph H
    # Note that rate = M/N
    # Design rate is 1 - j/k where j is the degree of variable nodes, k degree of check nodes
    # i.e. j ones per column, k ones per row. Rounding k up means no need to headache about running out of spots
    # Pick j = 3 from convention in the literature (smallest j with good theoretical properties)
    # So with M/N = 1 - 3/k, we pick k = ceil(3 / (1 - M/N))

    # However, we must have 3N ones in total, so we actually end up using a bigger k sometimes
    # Notably, we use the larger of k or 3N/M ones per row
    # Theoretically, this actually results in a better rate (though more expensive decoding process)

    # Alg is:
    # For each column, pick three random positions to put 1's s.t. no row ends up with more than k 1's in it
    # Inputs
    #  M        : num of rows (frequently k in literature, the dimension of the linear code)
    #  N        : num of columns
    #
    # Output
    #  H        : LDPC adjacency matrix   (N - k) x N = (N - M) x N

    j = 3
    M = N - M # dimension generated here is actually N - M
    H = falses(M, N)
    num_ones = zeros(Int64, M) #number of ones per row so far
    for current_col = 1:N
        # For uniformity, we want to have 3N/M ones in total at the end
        # So as we expand, say we've filled N' columns so far. We have 3N'
        # ones, so 3N'/M ones per row. So, to preserve uniformity as we go,
        # add to rows with less than ceil(3(N'+1)/M) ones where N' is the current
        # column we're filling. current_col is N' here.

        legal = [i for (i, val) in enumerate(num_ones) if val < ((3 * current_col / M))]
        if size(legal)[1] < j
            # Need j - size(legal) more spots, pick them uniformly randomly
            d = shuffle(collect(setdiff(Set(1:M), Set(legal))))
            additional_options = d[1:j - size(legal)[1]]
            append!(legal, additional_options)
        end

        choices = shuffle(legal)[1:j]
        for idx in choices
            H[idx, current_col] = 1
            num_ones[idx] += 1
        end
    end

    # An additional optimization could be to expand the girth, but expensive and probably unnecessary
    return H
end

function generate_LDPC_generator_matrix(H::BitArray)
    # https://en.wikipedia.org/wiki/Low-density_parity-check_code#Operational_use
    # has instructions for what we're doing, but essentially:
    # 1. Gaussian eliminate s.t. H = [P | I] as a block matrix
    # 2. Return G = [I | P^\top]
    # Returns nothing if H has no corresponding generator (rank-deficient)
    H = copy(H)
    M = size(H)[1]
    N = size(H)[2]
    pivot_row = 1
    pivot_col = N - M + 1
    

    for _ in 1:M
        # println("In pivot positon (", pivot_row, ", ", pivot_col, ")")
        # 1. Swap a row with a 1 into this pivot position if need be
        if !H[pivot_row, pivot_col]
            success = false
            for j in (pivot_row + 1):M
                if H[j, pivot_col]
                    # switch the found row and the pivot row
                    for k in 1:N
                        H[j, k], H[pivot_row, k] = H[pivot_row, k], H[j, k]
                    end
                    success = true
                    break
                end
            end
            if !success
                return nothing #Rank failure if we reach here
            end
        end
        
        # 2. Kill the 1s below this in this column
        for j in (pivot_row + 1):M
            if H[j, pivot_col]
                for k in 1:N
                    H[j, k] = xor(H[j, k], H[pivot_row, k])
                end
            end
        end

        # 3. Kill the 1s above this
        for j in 1:(pivot_row - 1)
            if H[j, pivot_col]
                for k in 1:N
                    H[j, k] = xor(H[j, k], H[pivot_row, k])
                end
            end
        end

        pivot_row += 1
        pivot_col += 1
    end

    # Extract P
    k = N - M
    P = transpose(H[:, 1:k])

    # Construct and return G
    G = falses(k, N)
    G[1:k, (k + 1):N] = P
    for i in 1:k
        G[i, i] = 1
    end
    return G
end

function generate_LDPC_code(M::Int64, N::Int64)
    H = nothing
    G = nothing
    # Sometimes random check matrices will have insufficient rank (no corresponding generator)
    # So we just resample until it's good. Ideally the probability of insufficient rank goes away for larger sizes,
    # something to mention/think about later than I am writing this.
    while G == nothing
        H = generate_LDPC_adjacency_matrix(M, N)
        G = generate_LDPC_generator_matrix(H)
    end
    return H, G
end

function check_pair(G, H)
    C = G * transpose(H)
    for i in 1:size(C)[1]
        for j in 1:size(C)[2]
            if C[i, j] % 2 != 0
                return false
            end
        end
    end
    return true
end

function do_test(M::Int64, N::Int64, iters::Int64)
    for _ in 1:iters
        H, G = generate_LDPC_code(M, N)
        println(size(H), size(G))
        if !check_pair(G, H)
            println("FAIL")
        end
    end
end

function encode_LDPC(input::BitVector, G::BitArray)
    return BitVector(vec((transpose(input) * G) .% 2))
end

function decode_LDPC(y::BitVector, H::BitArray, p::Float64, iters::Int64=30)
    # y is the received codeword of length N
    # H is the check matrix as usual of dimensions (N - k) x N
    N = size(H)[2]
    k = N - size(H)[1]
    M = N - k

    L = zeros(N) # LLR for the variable nodes
    lambda = zeros(M) # LLR for the check nodes
    messagesToCheckNodes = zeros(N, M) # Messages passed from variable nodes to check nodes. Omega matrix
    messagesToVarNodes = zeros(M, N) # Messages passed from check nodes to variable nodes. Lambda matrix

    # 1. Initialize the variable nodes based on received information (the identity part of the codeword)
    # We initialize each variable node's LLR
    # Assuming uniform prior, L_i = ln(p/(1-p)) when y_i = 0 and ln((1-p)/p) when y_i = 1
    lnq = log(p) - log(1-p)
    for i in 1:N
        L[i] = lnq * ((-1) ^ y[i])
    end
    println(L)

    # 30 rounds of decoding based on literature and practice
    for iter in 1:iters
        # 2. Update the check nodes and generate new check to variable messages
        for i in 1:M
            lambda[i] = 0
            for j in 1:N
                if H[i, j]
                    lambda[i] += log(tanh(abs(messagesToCheckNodes[j, i]) / 2))
                end
            end
        end

        println(lambda)
        for i in 1:M
            for j in 1:N
                Lambda[i, j] = 2 * atanh(exp(log(lambda[i]) - log(tanh(Omega[j, i] / 2))))
            end
        end

        # 3. Update the variable nodes and generate new variable to check messages
        for i in 1:N
            for j in 1:M
                if H[j, i]
                    L[i] += Lambda[j, i]
                end
            end
        end

        for i in 1:N
            for j in 1:M
                Omega[i, j] = L[i] - Lambda[j, i]
            end
        end
        
        # 4. Quantize and halt
        estimate = falses(N)
        for i in 1:N
            if L[i] >= 0
                estimate[i] = true
            end
        end

        if BitVector((transpose(estimate) * transpose(H)) .% 2) == falses(M)
            return estimate
        end
    end

    # Failure to decode
    return nothing
end

# H, G = generate_LDPC_code(15, 25)
# input = BitVector([1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0])
# codeword = encode_LDPC(input, G)
# recovered = decode_LDPC(codeword, H, 0.1)
# println(input)
# println(recovered)