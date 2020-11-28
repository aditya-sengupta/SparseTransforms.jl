using SparseTransforms
using Test
using Pipe.Pipe
using LinearAlgebra

function test_reconstruct()
    query_method = :simple
    delays_method = :nso

    n = 4
    locs = [4, 6, 10, 15]
    strengths = Float64.([2, 4, 1, 1])
    signal = TestSignal(n, locs, strengths, 0.0)
    b = get_b(signal; method=query_method)
    Ms = get_Ms(signal.n, b; method=query_method)
    D = get_D(signal.n; method=delays_method, num_delays=2*n)

    Us, Ss = Any[], Any[]
    K = binary_ints(signal.n)

    for M in Ms
        U, used_i = compute_delayed_subtransform(signal, M, D, fwht)
        U = hcat(U...)
        push!(Us, U)
        push!(Ss, (-1) .^(D * K))
    end

    select_froms = []
    for M in Ms
        selects = @pipe M' * K |> transpose |> Bool.(_) |> eachrow |> bin_to_dec.(_)
        push!(select_froms, selects)
    end

    cutoff = 2 * signal.noise_sd ^ 2 * (2 ^ (signal.n - b)) * (signal.n + 1) # noise threshold

    for method_name in [:nso] # all_methods["reconstruct"]
        println("Testing reconstruction with $method_name")
        singleton_guesses = []
        for (i, (U, S, select_from)) in enumerate(zip(Us, Ss, select_froms))
            col_gen = U |> eachrow |> enumerate
            for (j, col) in col_gen
                if col⋅col > cutoff
                    selection = findall(==(j-1), select_from) # pick all the k such that M.T @ k = j
                    k, sgn = singleton_detection(
                        col;
                        method=method_name,
                        selection=selection,
                        S_slice=S[:, selection],
                        n=signal.n
                    ) # find the best fit singleton
                    push!(singleton_guesses, k |> bin_to_dec)
                end
            end
        end
        @test singleton_guesses == [4, 6, 15, 6, 10, 15]
    end
end
