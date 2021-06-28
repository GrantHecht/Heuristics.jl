
abstract type Optimizer end

# Define general use solver 
function optimize!(opt::Optimizer, opts::Options)

    # Dispatch to _solve() with typeof opt
    res = _optimize!(opt, opts)

    #return results
    return res
end

function optimize!(opt::Optimizer)

    # Initialize options with defaults
    opts = Options()

    # Solve
    res = optimize!(opt, opts)

    # return results
    return res
end