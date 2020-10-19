"""
SPRIGHT decoding main file. Logic flow:

1. Generate a signal from input_signal.jl
2. Subsample from query.jl
3. Peel using reconstruct.jl
"""
module SPRIGHT

    export transform, method_test, method_report
    using ProgressMeter
    include("utils.jl")
    include("query.jl")
    include("reconstruct.jl")
    
    """
    Full SPRIGHT encoding and decoding. Implements Algorithms 1 and 2 from [2].
    (numbers) in the comments indicate equation numbers in [2].

    Arguments
    ---------
    signal : InputSignal object.
    The signal to be transformed / compared to.

    methods : Array{Symbol,1}
    The three symbols [query_method, delays_method, reconstruct_method].

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
    
    verbose : Bool
    Whether to print intermediate steps.

    report : Bool
    Whether to report on algorithm performance.

    Returns
    -------
    wht : Array{Float64,1}
    The WHT constructed by subsampling and peeling.
    """
    function transform(signal::InputSignal, methods::Array{Symbol,1}, verbose::Bool=false, report::Bool=False)
        query_method, delays_method, reconstruct_method = methods
        # check the condition for p_failure > eps
        # upper bound on the number of peeling rounds, error out after that point
        num_peeling = 0
        result = []
        wht = zeros(Float64, [signal.n])
        b = get_b(signal, method=query_method)
        Ms = get_Ms(signal.n, b, method=query_method)
        peeling_max = 2^b

        Us, Ss = [], []
        singletons = Dict()
        multitons = []
        if report
            used = Set()
        end
        if delays_method != :nso
            num_delays = signal.n + 1
        else
            num_delays = signal.n * Int64(log2(signal.n)) # idk
        end
        K = binary_ints(signal.n)

        # subsample, make the observation [U] and offset signature [S] matrices
        for M in Ms
            D = get_D(signal.n, method=delays_method; num_delays=num_delays)
            if verbose
                println("-----")
                println("a delay matrix")
                println(D)
            end
            U, used_i = compute_delayed_wht(signal, M, D)
            append!(Us, U)
            append!(Ss, (-1) .^(D ⋅ K))
            if report
                used = union(used, used_i)
            end
        end

        cutoff = 2 * signal.noise_sd ^ 2 * (2 ^ (signal.n - b)) * num_delays # noise threshold
        if verbose
            println("cutoff = ", cutoff)
        end

        # K is the binary representation of all integers from 0 to 2 ** n - 1.
        select_froms = map(M -> bin_to_dec.(M' ⋅ K)', Ms)
        # `select_froms` is the collection of 'j' values and associated indices
        # so that we can quickly choose from the coefficient locations such that M.T @ k = j as in (20)
        # example: ball j goes to bin at "select_froms[i][j]"" in stage i
        
        # begin peeling
        # index convention for peeling: 'i' goes over all M/U/S values
        # i.e. it refers to the index of the subsampling group (zero-indexed - off by one from the paper).
        # 'j' goes over all columns of the WHT subsample matrix, going from 0 to 2 ** b - 1.
        # e.g. (i, j) = (0, 2) refers to subsampling group 0, and aliased bin 2 (10 in binary)
        # which in the example of section 3.2 is the multiton X[0110] + X[1010] + W1[10]
        
        # a multiton will just store the (i, j)s in a list
        # a singleton will map from the (i, j)s to the true (binary) values k.
        # e.g. the singleton (0, 0), which in the example of section 3.2 is X[0100] + W1[00]
        # would be stored as the dictionary entry (0, 0): array([0, 1, 0, 0]).
        
        multitons_found = True
        iters = 0
        max_iters = 2 ^ (b + 1)
        while multitons_found && (num_peeling < peeling_max) && (iters < max_iters)
            if verbose
                println("-----")
                println("the measurement matrix")
                for U in Us
                    println(U)
                end
            end
            
            # first step: find all the singletons and multitons.
            singletons = Dict() # dictionary from (i, j) values to the true index of the singleton, k.
            multitons = [] # list of (i, j) values indicating where multitons are.
            
            for (i, (U, S, select_from)) in enumerate(zip(Us, Ss, select_froms))
                for (j, col) in enumerate(U.T)
                    if col⋅col > cutoff
                        selection = findall(==(j), select_from) # pick all the k such that M.T @ k = j
                        k, sgn = singleton_detection(
                            col, 
                            method=self.reconstruct_method; 
                            selection=selection, 
                            S_slice=S[:, selection], 
                            n=signal.n
                        ) # find the best fit singleton
                        k_dec = bin_to_dec(k)
                        rho = (S[:,k_dec] ⋅ col) * sgn / length(col)                    
                        residual = col - sgn * rho * S[:,k_dec] 
                        if verbose
                            println((i, j), np.inner(residual, residual))
                        end
                        if residual ⋅ residual > cutoff
                            append!(multitons, (i, j))
                        else # declare as singleton
                            singletons[(i, j)] = (k, rho, sgn)
                            if verbose
                                print("amplitude, ", rho)
                            end
                        end # if residual norm > cutoff
                    end # if col norm > cutoff
                end # for col
            end # for U, S, select_from
                            
            # all singletons and multitons are discovered
            if verbose
                println("singletons:")
                for (ston_key, ston_value) in singletons
                    println(ston_key, bin_to_dec(ston_value[1]))
                end

                println("Multitons : ", multitons)
            end
            
            # WARNING: this is not a correct thing to do
            # in the last iteration of peeling, everything will be singletons and there
            # will be no multitons
            if length(multitons) == 0 # no more multitons, and can construct final WHT
                multitons_found = false
            end
                
            # balls to peel
            balls_to_peel = Set()
            ball_values = Dict()
            ball_sgn = Dict()
            for (_, (k, rho, sgn)) in singletons
                push!(balls_to_peel, bin_to_dec(k))
                ball_values[ball] = rho
                ball_sgn[ball] = sgn
            end
                
            if verbose
                println("these balls will be peeled")
                print(balls_to_peel)
            end
            # peel
            iters += 1
            for ball in balls_to_peel
                num_peeling += 1
                k = dec_to_bin(ball, signal.n)
                potential_peels = [(l, bin_to_dec(M.T.dot(k))) for (l, M) in enumerate(Ms)]

                append!(result, (k, ball_sgn[ball]*ball_values[ball]))
                
                for peel in potential_peels
                    signature_in_stage = Ss[peel[0]][:,ball]
                    to_subtract = ball_sgn[ball] * ball_values[ball] * signature_in_stage
                    Us[peel[0]][:,peel[1]] -= to_subtract
                    if verbose
                        println("this is subtracted:")
                        println(to_subtract)
                        println("Peeled ball ", bin_to_dec(k), " off bin ", peel)
                    end
                end # for peel
            end # for ball
        end # while
            
        loc = Set()
        for (k, value) in result # iterating over (i, j)s
            idx = bin_to_dec(k) # converting 'k's of singletons to decimals
            push!(loc, idx)
            if wht[idx] == 0
                wht[idx] = value
            else
                wht[idx] = (wht[idx] + value) / 2 
                # average out noise; e.g. in the example in 3.2, U1[11] and U2[11] are the same singleton,
                # so averaging them reduces the effect of noise.
            end
        end
        
        wht /= 2 ^ (signal.n - b)
        if report
            return wht, len(used), loc
        else
            return wht
        end
    end

    """
    Tests a method on a signal and reports its average execution time and sample efficiency.
    """
    function method_test(signal, num_runs=10)
        time_start = time()
        samples = 0
        successes = 0
        for i in ProgressBar(num_runs)
            wht, num_samples, loc = transform(signal, report=True)
            if loc == set(signal.loc)
                successes += 1
            end
            samples += num_samples
        end
        return (time() - time_start) / num_runs, successes / num_runs, samples / (num_runs * 2 ^ signal.n)
    end

    """
    Reports the results of a method_test.
    """
    function method_report(signal::InputSignal, num_runs::Int64=10)
        println("Testing SPRIGHT with query method $query_method, delays method $delays_method, reconstruct method $reconstruct_method.")
        t, s, sam = method_test(signal, num_runs)
        print("Average time in seconds: ", format(t))
        print("Success ratio: ", s)
        print("Average sample ratio: ", sam)
    end

end # module
