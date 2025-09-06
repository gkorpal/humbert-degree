#include("./expt/getDeg.jl")
#using Pkg
#Pkg.develop(path="../Humbert")
using Oscar # is_probable_prime
using Humbert

function getDeg(N::Int, n::Int, m::Int)
    primes = primeList(N, n, m)
    last_p = primes[end]
    output_filename = "deg_$(p)_$(last_p).txt"
    open(output_filename, "w") do io
        for p in primes
            println(io, "p = ", p)
            ell = getEll(p, 1)
            println(io, "ℓ = ", ell)
            try
                NDict = minDeg(p, ell, 2*ell*p, 2*ell*p+p)
                println(io, "max(min deg) = ", maximum(keys(NDict)))
                println(io, "minimum degree frequency distribution")
                for (n, freq) in sort(collect(NDict), by = x -> x[2], rev = true)
                    println(io, n, " => ", freq)
                end
                println(io, "")   
            catch e 
                println(io, "skipping $p since file doesn't exist.\n")
            end
        end
    end
    return nothing
end

