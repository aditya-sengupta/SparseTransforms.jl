"""
Sparse transform (FFAST/SPRIGHT) decoding main file. Logic flow:

1. Generate a signal from input_signal.jl
2. Subsample from query.jl
3. Peel using reconstruct.jl
"""
module SparseTransforms
    export all_methods, spright, ffast, method_test, method_report
    using ProgressMeter
    include("reconstruct.jl")
    export fwht, bin_to_dec, dec_to_bin, binary_ints, sign_spright, expected_bin
    export Signal, TestSignal, InputSignal
    export get_D, get_b, get_Ms, subsample_indices, compute_delayed_subtransform
    export singleton_detection, bin_cardinality

    all_methods = Dict(
        "query" => [:simple],
        "delays" => [:identity_like, :random, :nso],
        "reconstruct" => [:noiseless, :mle, :nso]
    )

    function spright(signal::Signal, methods::Array{Symbol,1}; verbose::Bool=false, report::Bool=false)
        return transform(signal, methods, fwht; verbose=verbose, report=report)
    end

    function ffast(signal::Signal, methods::Array{Symbol,1}; verbose::Bool=false, report::Bool=false)
        return transform(signal, methods, fft; verbose=verbose, report=report)
    end

    """
    Full SPRIGHT encoding and decoding. Implements Algorithms 1 and 2 from [2].
    (numbers) in the comments indicate equation numbers in [2].

    Arguments
    ---------
    signal : TestSignal object.
    The signal to be transformed / compared to.

    methods : Array{Symbol,1}
    The three symbols [query_method, delays_method, reconstruct_method].
    All implemented methods are available in `all_methods`: if you add a new method, make sure to update `all_methods`.

        query_method : Symbol
        The method to generate the sparsity coefficient and the subsampling matrices.
        Currently implemented methods:
            :simple : choose some predetermined matrices based on problem size.

        delays_method : Symbol
        The method to generate the matrix of delays.
        Currently implemented methods:
            :identity_like : return a zero row and the identity matrix vertically stacked together.
            :random : make a random set of delays.
            :nso : set delays according to the NSO-SPRIGHT algorithm.

        reconstruct_method : Symbol
        The method to detect singletons.
        Currently implemented methods:
            :noiseless : decode according to [2], section 4.2, with the assumption the signal is noiseless.
            :mle : naive noisy decoding; decode by taking the maximum-likelihood singleton that could be at that bin.
            :nso : reconstruct according to the NSO-SPRIGHT algorithm.

    transform : Function
    The base transform: either `fwht` or `fft`.

    verbose : Bool
    Whether to print intermediate steps.

    report : Bool
    Whether to report on algorithm performance.

    Returns
    -------
    wht : Array{Float64,1}
    The WHT constructed by subsampling and peeling.
    """
    function transform(signal::Signal, methods::Array{Symbol,1}, transform::Function; verbose::Bool=false, report::Bool=false)
        for (method_type, method_name) in zip(["query", "delays", "reconstruct"], methods)
            impl_methods = all_methods[method_type]
            @assert method_name in impl_methods "$method_type method must be one of $impl_methods"
        end
        query_method, delays_method, reconstruct_method = methods
        # check the condition for p_failure > eps
        # upper bound on the number of peeling rounds, error out after that point

        num_peeling = 0
        locs = []
        strengths = []
        wht = Dict()
        b = get_b(signal; method=query_method)
        Ms = get_Ms(signal.n, b; method=query_method)
        peeling_max = 2^b

        Us = []
        if report
            used = Set()
        end

        if delays_method != :nso
            num_delays = signal.n + 1
        else
            num_delays = signal.n * Int64(ceil(log2(signal.n))) # idk
        end
        D = get_D(signal.n; method=delays_method, num_delays=num_delays)
        if reconstruct_method == :mle
            K = binary_ints(signal.n)
            S = (-1) .^(D * K)
        end

        # subsample, make the observation [U] matrices
        for M in Ms
            U, used_i = compute_delayed_subtransform(signal, M, D, transform)
            push!(Us, U)
            if report
                used = union(used, used_i)
            end
        end

        cutoff = 4 * signal.noise_sd ^ 2 * (2 ^ (signal.n - b)) * num_delays # noise threshold

        # K is the binary representation of all integers from 0 to 2 ** n - 1.
        select_froms = []
        for M in Ms
            selects = 0
            if reconstruct_method == :mle
                selects = @pipe M' * K |> transpose |> Bool.(_) |> eachrow |> bin_to_dec.(_)
            end
            push!(select_froms, selects)
        end
        # `select_froms` is the collection of 'j' values and associated indices
        # so that we can quickly choose from the coefficient locations such that M.T * k = j as in (20)
        # example: ball j goes to bin at "select_froms[i][j]"" in stage i

        # begin peeling
        # index convention for peeling: 'i' goes over all M/U values
        # i.e. it refers to the index of the subsampling group (zero-indexed - off by one from the paper).
        # 'j' goes over all columns of the WHT subsample matrix, going from 0 to 2 ** b - 1.
        # e.g. (i, j) = (0, 2) refers to subsampling group 0, and aliased bin 2 (10 in binary)
        # which in the example of section 3.2 is the multiton X[0110] + X[1010] + W1[10]

        # a multiton will just store the (i, j)s in a list
        # a singleton will map from the (i, j)s to the true (binary) values k.
        # e.g. the singleton (0, 0), which in the example of section 3.2 is X[0100] + W1[00]
        # would be stored as the dictionary entry (0, 0): array([0, 1, 0, 0]).

        multitons_found = true
        iters = 0
        max_iters = 2 ^ (b + 1)
        while multitons_found && (num_peeling < peeling_max) && (iters < max_iters)

            # first step: find all the singletons and multitons.
            singletons = Dict() # dictionary from (i, j) values to the true index of the singleton, k.
            multitons = [] # list of (i, j) values indicating where multitons are.

            for (i, (M, U, select_from)) in enumerate(zip(Ms, Us, select_froms))
                col_gen = U |> eachrow |> enumerate
                for (j, col) in col_gen
                    if col⋅col > cutoff
                        selection = findall(==(j), select_from) # pick all the k such that M.T @ k = j
                        slice = 0
                        if reconstruct_method == :mle
                            slice = S[:, selection]
                        end
                        k, sgn = singleton_detection(
                            col;
                            method=reconstruct_method,
                            selection=selection,
                            S_slice=slice,
                            n=signal.n
                        ) # find the best fit singleton
                        s_k = (-1) .^ (D * k)
                        ρ = (s_k ⋅ col) * sgn / length(col)
                        residual = col - sgn * ρ * s_k
                        if (expected_bin(k, M) != j) || (residual ⋅ residual > cutoff)
                            push!(multitons, [i, j])
                        else # declare as singleton
                            singletons[(i, j)] = (k, ρ, sgn)
                        end # if residual norm > cutoff
                    end # if col norm > cutoff
                end # for col
            end # for U, select_from

            # all singletons and multitons are discovered
            if verbose
                println("Singletons:")
                for (ston_key, ston_value) in singletons
                    println(ston_key, bin_to_dec(ston_value[1]))
                end

                println("Multitons : $multitons")
            end

            # WARNING: this is not a correct thing to do
            # in the last iteration of peeling, everything will be singletons and there
            # will be no multitons
            if length(multitons) == 0 # no more multitons, and can construct final transform
                multitons_found = false
            end

            # balls to peel
            balls_to_peel = Set()
            ball_values = Dict()
            ball_sgn = Dict()
            for (_, (k, rho, sgn)) in singletons
                ball = bin_to_dec(k)
                push!(balls_to_peel, ball)
                ball_values[ball] = rho
                ball_sgn[ball] = sgn
            end

            if verbose
                println("Balls to be peeled")
                println(balls_to_peel)
            end
            # peel
            iters += 1
            for ball in balls_to_peel
                num_peeling += 1
                k = dec_to_bin(ball, signal.n)

                push!(locs, k)
                push!(strengths, ball_sgn[ball]*ball_values[ball])

                for (l, M) in enumerate(Ms)
                    peel = @pipe M' * k |> Bool.(_) |> bin_to_dec |> _ + 1
                    signature_in_stage = (-1) .^ (D * k)
                    to_subtract = ball_sgn[ball] * ball_values[ball] * signature_in_stage
                    U = Us[l]
                    U_slice = U[peel, :]
                    Us[l][peel,:] = U_slice - to_subtract
                    if verbose
                        println("Peeled ball $(bin_to_dec(k)) off bin ($l, $peel)")
                    end
                end # for peel
            end # for ball
        end # while

        loc = Set()
        for (k, value) in zip(locs, strengths) # iterating over (i, j)s
            idx = bin_to_dec(k) # converting 'k's of singletons to decimals
            push!(loc, idx)
            if !haskey(wht, idx+1)
                wht[idx] = value
            end
        end

        wht = Dict(i => x / (2 ^ b) for (i,x) in wht)
        if report
            return wht, length(used)
        else
            return wht
        end
    end

    """
    Tests a method on a signal and reports its average execution time and sample efficiency.
    """
    function method_test(signal::TestSignal, methods::Array{Symbol,1}; num_runs::Int64=10)
        time_start = time()
        num_samples = 0
        num_successes = 0
        for i in ProgressBar(num_runs)
            wht, samples = spright(signal, methods; report=True)
            success = length(signal.loc) == length(spright_wht)
            for (l, s) in zip(signal.loc, signal.strengths)
                success = success & isapprox(wht[l], s, atol=5*signal.σ)
            end
            if success
                num_successes += 1
            end
            num_samples += samples
        end
        return (time() - time_start) / num_runs, num_successes / num_runs, num_samples / (num_runs * 2 ^ signal.n)
    end

    """
    Reports the results of a method_test.
    """
    function method_report(signal::TestSignal, num_runs::Int64=10)
        println("Testing SPRIGHT with query method $query_method, delays method $delays_method, reconstruct method $reconstruct_method.")
        t, s, sam = method_test(signal, num_runs)
        print("Average time in seconds: ", format(t))
        print("Success ratio: ", s)
        print("Average sample ratio: ", sam)
    end

end # module
