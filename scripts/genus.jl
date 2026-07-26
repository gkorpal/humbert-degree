using Oscar

# q_5 = <1> ⊥ 4(b_0(t1,t2) ⊥ b_0(t3,t4)), b_0 = x^2+xy+(p+1)/4 y^2 the norm
# form of Z[(1+√-p)/2], p ≡ 3 (mod 4). Splitting off the norm-1 vector and
# then descaling by 4 gives a bijection
#     gen(q)  --(×4)-->  gen(4q)  --(<1> ⊥ ·)-->  {L ∈ gen(q_5) : min(L) = 1},
# so only the rank-4 form q = b_0 ⊥ b_0 need be enumerated. Its Gram matrix
# is genuinely half-integral (b_0 has odd middle coefficient), so it must be
# built over QQ directly rather than by halving a doubled integer matrix.

"""
    poly_form(M::ZZMatrix)

Quadratic form `1/2 * x^T M x` as a polynomial, for a doubled Gram matrix
`M` (M[i,i] = 2Q(e_i), M[i,j] = 2B(e_i,e_j) for i<j).
"""
function poly_form(M::ZZMatrix)
    n = number_of_rows(M)
    R, x = polynomial_ring(ZZ, n)
    f = R(0)
    for i in 1:n
        f += divexact(M[i, i], 2) * x[i]^2
        for j in (i+1):n
            f += M[i, j] * x[i] * x[j]
        end
    end
    return f
end

"""
    binary_gram(a::Integer, b::Integer, c::Integer) -> QQMatrix

Gram matrix `[a b/2; b/2 c]` of `a*x^2 + b*x*y + c*y^2`, exact for odd `b`.
"""
binary_gram(a::Integer, b::Integer, c::Integer) =
    matrix(QQ, 2, 2, [QQ(a), QQ(b, 2), QQ(b, 2), QQ(c)])

"""
    quaternary_form_q(p::Integer) -> QQMatrix

Gram matrix of `q := b_0 ⊥ b_0`, the principal binary form of discriminant
`-p` taken twice.
"""
function quaternary_form_q(p::Integer)
    G = binary_gram(1, 1, div(p + 1, 4))
    return block_diagonal_matrix([G, G])
end

"""
    rational_to_integer_matrix(M::QQMatrix) -> ZZMatrix

Convert an integral `QQMatrix` to a `ZZMatrix`, asserting integrality
rather than truncating.
"""
function rational_to_integer_matrix(M::QQMatrix)
    n = number_of_rows(M)
    m = number_of_columns(M)
    ZM = zero_matrix(ZZ, n, m)
    for i in 1:n, j in 1:m
        x = M[i, j]
        @assert isone(denominator(x)) "entry ($i,$j) = $x is not an integer"
        ZM[i, j] = numerator(x)
    end
    return ZM
end

"""
    genus5(p::Integer)

Compute representatives of the genus of `q_5` with minimum 1, for a prime
`p ≡ 11 (mod 12)`, and write them as quinary polynomials to `Gen5_\$p.txt`.
"""
function genus5(p::Integer)
    println("p = $p")

    Gq = quaternary_form_q(p)
    Lq = integer_lattice(; gram = Gq)
    @assert gram_matrix(Lq) == Gq "Lq's Gram matrix does not match q"

    reps = genus_representatives(Lq)

    output_filename = "Gen5_$(p).txt"
    open(output_filename, "w") do io
        println(io, "p = ", p)
        println(io, "method: genus of q = b_0 ⊥ b_0 (rank 4, disc -p each block), scaled by 4\n")
        for L in reps
            M4 = rational_to_integer_matrix(4 * gram_matrix(L))
            Mfull = block_diagonal_matrix([ZZ[2;;], M4 .* 2])
            println(io, "q = ", poly_form(Mfull), "\n")
        end
        println(io, "Total Gen5 forms rep 1: ", length(reps))
    end
    println("  -> wrote $(length(reps)) forms to $output_filename\n")
    return nothing
end

# for p in 11:12:660
#     is_probable_prime(p) && genus5(p)
# end

# Usage: julia genus.jl 503
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || (println("Usage: julia genus.jl <prime>"); exit(1))

    p = parse(Int, ARGS[1])
    p % 12 == 11 && is_probable_prime(p) ||
        error("Expected a prime congruent to 11 mod 12.")

    genus5(p)
end