#using Pkg
#Pkg.develop(path="../Humbert")
using Humbert

function getRHIs(N::Int, n::Int, m::Int)
    primes = primeList(N, n, m)
    for p in primes
        ell = getEll(p, 1)
        a = 2 * ell * p
        b = a + p
        allRHI(p, ell, a, b)
    end
    return nothing
end

getRHIs(19,1,3) # force compile
