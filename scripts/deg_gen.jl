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

const DEGREE2_MONOMIALS = [
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
    parse_gen_file(filename)

Parse a `Gen5_p.txt` file: extract the prime `p` and, for each `q = ...`
line, recover the degree form `(q - x1^2) / 4`. Blank lines, a leading
`method: ...` comment, and a trailing `Total Gen5 forms rep N: count` line
are recognized and skipped; the trailing count, if present, is checked
against the number of forms parsed. Returns `(p, forms)` where `forms` maps
type index to `(q, deg)`.
"""
function parse_gen_file(filename::String)
    prime_re = r"^\s*p\s*=\s*(\d+)\s*$"
    q_re = r"^\s*q\s*=\s*(.+?)\s*$"
    total_re = r"^\s*Total\s+Gen5\s+forms\s+rep\s+\d+:\s*(\d+)\s*$"
    method_re = r"^\s*method\s*:"

    prime = nothing
    expected_total = nothing
    forms = Dict{Int,NTuple{2,PolyElem}}()
    type_index = 0

    for line in eachline(filename)
        isempty(strip(line)) && continue

        if (m = match(prime_re, line)) !== nothing
            prime = parse(Int, m.captures[1])
        elseif (m = match(q_re, line)) !== nothing
            q = parse_poly(m.captures[1])
            deg = divexact(q - x1^2, 4)
            type_index += 1
            forms[type_index] = (q, deg)
        elseif (m = match(total_re, line)) !== nothing
            expected_total = parse(Int, m.captures[1])
        elseif match(method_re, line) === nothing
            @warn "unrecognized line in $filename" line
        end
    end

    if expected_total !== nothing
        @assert length(forms) == expected_total "parsed $(length(forms)) forms, file declares $expected_total"
    end

    return prime, forms
end

"""
    min_deg_distribution(p) -> Dict{Int,Int}

Read `Gen5_p.txt`, compute the minimum vector of the lattice for each
degree form's Gram matrix, and return the frequency distribution of minima.
"""
function min_deg_distribution(p::Integer)
    prime, forms = parse_gen_file("Gen5_$p.txt")
    @assert prime == p "file mismatch"

    freq = Dict{Int,Int}()
    for (_, deg) in values(forms)
        c = [coeff(deg, mono) for mono in DEGREE2_MONOMIALS]
        A = ZZ[2c[1] c[2] c[3] c[4]; c[2] 2c[5] c[6] c[7]; c[3] c[6] 2c[8] c[9]; c[4] c[7] c[9] 2c[10]]
        n = Int(minimum(integer_lattice(gram = A .÷ 2)))
        freq[n] = get(freq, n, 0) + 1
    end
    return freq
end

"""
    get_deg(p)

Write the minimum-degree frequency distribution for prime `p` to `deg_p.txt`.
If `Gen5_p.txt` is missing or parsing fails, print a message instead of
creating an output file.
"""
function get_deg(p::Integer)
    freq = try
        min_deg_distribution(p)
    catch e
        if e isa SystemError || e isa Base.IOError
            println("skipping $p: file doesn't exist")
        else
            println("skipping $p: ", sprint(showerror, e))
        end
        return nothing
    end

    open("deg_gen_$p.txt", "w") do io
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

for p in 11:12:660
    is_probable_prime(p) && get_deg(p)
end