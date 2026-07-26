using Oscar

const R, (x1, x2, x3, x4, x5) = polynomial_ring(ZZ, 5)
const PolyElem = elem_type(R)
const VARS = Dict(:x1 => x1, :x2 => x2, :x3 => x3, :x4 => x4, :x5 => x5)

"""
    poly_from_expr(ex)

Evaluate a parsed arithmetic expression (`+`, `-`, `*`, `^`, integers, and
`x1`..`x5`) into an element of `R`, without calling `eval`.
"""
poly_from_expr(ex::Integer) = R(ex)
poly_from_expr(ex::Symbol) = VARS[ex]
function poly_from_expr(ex::Expr)
    ex.head == :call || error("unsupported expression: $ex")
    op = ex.args[1]
    op === :^ && return poly_from_expr(ex.args[2])^(ex.args[3]::Integer)
    args = poly_from_expr.(ex.args[2:end])
    op === :+ && return reduce(+, args)
    op === :- && return length(args) == 1 ? -args[1] : reduce(-, args)
    op === :* && return reduce(*, args)
    error("unsupported operator: $op")
end

"""
    parse_poly(str)

Parse a polynomial expression string in `x1`..`x5` into an element of `R`.
"""
parse_poly(str::AbstractString) = poly_from_expr(Meta.parse(str))

"""
    parse_file(filename) -> (prime, types_data)

Parse `filename`, returning the prime `p` and a dict mapping each type number
to its `(q, deg)` pair.
"""
function parse_file(filename::AbstractString)
    prime = nothing
    types_data = Dict{Int,NTuple{2,PolyElem}}()

    current_type = 0
    current_q = nothing
    current_deg = nothing

    prime_re = r"^\s*p\s*=\s*(\d+)"
    type_re = r"^\s*Type\s+(\d+)"
    q_re = r"^\s*q\(ExE,θ\)\s*=\s*(.+)"
    deg_re = r"^\s*deg\(ExE,θ\)\s*=\s*(.+)"

    for line in eachline(filename)
        m = match(prime_re, line)
        if m !== nothing
            prime = parse(Int, m.captures[1])
            continue
        end

        m = match(type_re, line)
        if m !== nothing
            if current_type != 0
                types_data[current_type] = (current_q, current_deg)
            end
            current_type = parse(Int, m.captures[1])
            current_q = nothing
            current_deg = nothing
            continue
        end

        if current_type != 0
            m = match(q_re, line)
            if m !== nothing
                current_q = parse_poly(strip(m.captures[1]))
                continue
            end

            m = match(deg_re, line)
            if m !== nothing
                current_deg = parse_poly(strip(m.captures[1]))
                continue
            end
        end
    end

    if current_type != 0
        types_data[current_type] = (current_q, current_deg)
    end

    return prime, types_data
end

const DEGREE2_MONOMIALS = [
    [2, 0, 0, 0, 0],
    [1, 1, 0, 0, 0],
    [1, 0, 1, 0, 0],
    [1, 0, 0, 1, 0],
    [1, 0, 0, 0, 1],
    [0, 2, 0, 0, 0],
    [0, 1, 1, 0, 0],
    [0, 1, 0, 1, 0],
    [0, 1, 0, 0, 1],
    [0, 0, 2, 0, 0],
    [0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1],
    [0, 0, 0, 2, 0],
    [0, 0, 0, 1, 1],
    [0, 0, 0, 0, 2],
]

"""
    min_deg(p) -> Dict{Int,Int}

Read the degree forms from "RHI1_p.txt", compute each lattice's minimum
vector from its gram matrix, and return the frequency distribution of those
minima.
"""
function min_deg(p::Integer)
    filename = "RHI1_$(p).txt"
    prime, forms = parse_file(filename)
    @assert prime == p "file mismatch"

    freq = Dict{Int,Int}()

    for (_, q) in values(forms)
        c = [coeff(q, mono) for mono in DEGREE2_MONOMIALS]

        gram = ZZ[
            2*c[1] c[2] c[3] c[4]
            c[2] 2*c[6] c[7] c[8]
            c[3] c[7] 2*c[10] c[11]
            c[4] c[8] c[11] 2*c[13]
        ]

        L = integer_lattice(gram = gram .÷ 2)
        n = Int(minimum(L))
        freq[n] = get(freq, n, 0) + 1
    end
    return freq
end

"""
    get_deg(p)

Write `min_deg(p)`'s minimum-vector frequency distribution to `deg_p.txt`.
If `RHI1_p.txt` is missing or parsing fails, print a message instead of
creating an output file.
"""
function get_deg(p::Integer)
    freq = try
        min_deg(p)
    catch e
        if e isa SystemError || e isa Base.IOError
            println("skipping $p: file doesn't exist")
        else
            println("skipping $p: ", sprint(showerror, e))
        end
        return nothing
    end

    open("deg_$(p).txt", "w") do io
        println(io, "p = ", p)
        println(io, "max(min deg) = ", maximum(keys(freq)))
        println(io, "minimum degree frequency distribution")
        for (n, count) in sort(collect(freq), by = last, rev = true)
            println(io, n, " => ", count)
        end
        println(io)
    end
    return nothing
end

# for p in 11:12:660
#     is_probable_prime(p) && get_deg(p)
# end