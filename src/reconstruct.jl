"""
Methods for the reconstruction engine; specifically, to

1. carry out singleton detection
2. get the cardinalities of all bins in a subsampling group (debugging only).
"""

include("query.jl")

"""
MLE value estimation of a singleton given its frequency.
"""
function estimate_value(k::BitArray{1}, D::BitArray{2}, U_slice)
    s_k = (-1) .^ (D * k)
    return (s_k ⋅ U_slice) / length(U_slice), s_k
end


"""
Noiseless-case singleton detection. Assumes P = n + 1 and D = [0; I].
See singleton_detection for full signature.
"""
function singleton_detection_noiseless(U_slice; kwargs...)
    s_k = U_slice[2:length(U_slice)] / U_slice[1]
    return -sign.(s_k) .== 1, U_slice[1], s_k
end

"""
Singleton detection that uses MLE: looks at the residuals created by peeling off each possible singleton.
See singleton_detection for full signature.
"""
function singleton_detection_mle(U_slice; kwargs...)
    selection, S_slice, n = kwargs[:selection], kwargs[:S_slice], kwargs[:n]
    P = size(S_slice, 1)
    alphas = (1 / P) * S_slice' * U_slice
    residual_vecs = (S_slice * Diagonal(alphas)) .- U_slice
    residuals = dropdims(sum(abs2, residual_vecs, dims=1); dims=1)
    k_sel = argmin(residuals)
    return dec_to_bin(selection[k_sel] - 1, n), alphas[k_sel], S_slice[k_sel]
end

"""
Non-sample-optimal singleton detection.
"""
function singleton_detection_nso(U_slice; kwargs...)
    n = kwargs[:n]
    D = kwargs[:D]
    p1 = length(U_slice) ÷ (n + 1)
    chunks = @pipe reshape(U_slice, (n + 1, p1)) |> transpose |> sign_spright.(_)
    chunks = chunks[:,2:n+1] .⊻ chunks[:,1]
    est_k = @pipe sum(chunks, dims=1) |> dropdims(_, dims=1) |> (_ .> (p1 ÷ 2))
    return est_k, estimate_value(est_k, D, U_slice)...
end

"""
Robust sign estimation for sample-optimal singleton detection.
"""
function estimate_sign(U_slice)
    return sum(sign_spright.(U_slice)) > length(U_slice) ÷ 2
end

"""
Sample-optimal singleton detection.
"""
function singleton_detection_so(U_slice; kwargs...)
    n = kwargs[:n]
    p1, p2, p3 = kwargs[:num_delays]
    code = kwargs[:code]
    sgn = estimate_sign(U_slice[:, p1+1:p1+p2])
    k = decode_with(U_slice[:, p1+p2+1:p1+p2+p3], sgn, code)
    return k, sgn
end

singleton_detection_lookup = Dict(
    :mle => singleton_detection_mle,
    :noiseless => singleton_detection_noiseless,
    :nso => singleton_detection_nso,
    :so => singleton_detection_so
)

"""
Finds the true index of a singleton, or the best-approximation singleton of a multiton.

Arguments
---------
U_slice: P-element Array{Float64,1}
The WHT component of a subsampled bin, with element i corresponding to delay i.

Returns
-------
k : n-element Array{Bool,1}
Index of the corresponding right node.
"""
function singleton_detection(U_slice; method=:mle, kwargs...)
    return singleton_detection_lookup[method](U_slice; kwargs...)
end

"""
Computes delayed WHT observations and declares cardinality based on that.
2 is a stand-in for any cardinality > 1. For debugging purposes.

Arguments
---------
signal : TestSignal
The input signal object.

M : n×b Array{Bool,2}
The subsampling matrix; takes on binary values.

D : num_delays×n Array{Int64,2}
The delays matrix.

Returns
-------
(these are still as they are in Python - need to check dims first)
cardinality : numpy.ndarray
0 or 1 if the bin is a zeroton or singleton resp.; 2 if multiton.

singleton_indices : list
A list (in decimal form for compactness) of the k values of the singletons.
Length matches the number of 1s in cardinality.
"""
function bin_cardinality(signal::Signal, M, D; method=:mle, verbose=false)
    b = size(M, 2)
    U, inds = compute_delayed_wht(signal, M, D)
    cardinality = ones(Int64, signal.n) # vector of indicators
    singleton_indices = []
    cutoff = 2 * signal.noise_sd^2 * (2^(signal.n - b)) * size(D, 1)
    if signal.noise_sd > 0
        K = binary_ints(signal.n)
        S = (-1).^ (D * K)
    end
    col_generator = U |> eachrow |> enumerate
    for (i, col) in col_generator
        sgn = 1
        if verbose
            println("Column:   ", col)
        end
        if col ⋅ col <= cutoff
            cardinality[i] = 0
        else
            if signal.noise_sd == 0
                k, sgn = singleton_detection_noiseless(col)
            else
                selection = @pipe bin_to_dec.(M' * k)' |> findall(==(i), _)
                k, sgn = singleton_detection(col, method=method, selection=selection, S_slice=S[:, selection], n=signal.n)
            end
            rho = @pipe col |> abs.(_) |> mean
            residual = col - sgn * rho * (-1).^ (D * k)
            if verbose
                println("Residual: ", residual)
            end
            if residual ⋅ residual > cutoff
                cardinality[i] = 2
            else
                append!(singleton_indices, bin_to_dec(k))
            end
        end
    end
    return cardinality, singleton_indices
end
