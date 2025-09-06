using Humbert

function getHh(N::Int, n::Int, m::Int)
    output_filename = "class_numbers_$(N)_$(n)_$(m).txt"
    primes = primeList(N, n, m)
    open(output_filename, "w") do io
        for p in primes
            println(io, "p = ", p)
            println(io, "[h, h(h+1)/2, H] = ", classNumbers(p), "\n")
        end
    end
    return nothing
end
