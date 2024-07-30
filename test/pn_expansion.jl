@testitem "PNTerm algebra" begin
    using DoubleFloats
    using PostNewtonian: PNExpansion

    for T ∈ [Float64, Float16, Double64]
        pn = rand(BBH; v=T(1) / 5, PNOrder=9//2)
        c = PNExpansionParameter(pn)
        z = zero(T)
        x = T(6) / 5
        y = T(17) / 5
        w = T(28) / 5
        v = T(27) / 13

        # Test behavior of `c` as the basic PNTerm
        for (term, c⁻¹exponent, coeff) ∈ (
            (c, -1, 1),
            (c^2, -2, 1),
            (x * c^2, -2, x),
            (c^2 * x, -2, x),
            (c^2 / x, -2, 1 / x),
            ((x * c)^2, -2, x^2),
            ((x * c^2) / c^4, 2, x),
            ((x * c^2) / c^-4, -6, x),
            ((x / c^2) / c^4, 6, x),
            ((x / c^2) * (y / c^4), 6, x * y),
        )
            @test PostNewtonian.c⁻¹exp(term) == c⁻¹exponent
            @test term.coeff == coeff
            @test term.coeff isa eltype(pn)
            @test length(term) == 1
            @test eltype(term) === T
        end

        # Test PNExpressions
        @test_throws ArgumentError PNExpansion((), 0)
        for (expr, expected) ∈ (
            (w - (x - y / c), (w - x, y)),
            (w - (x / c - y), (w + y, -x)),
            (x - y / c, (x, -y)),
            (x / c - y, (-y, x)),
            (x / c^6 + y / c, (z, y, z, z, z, z, x)),
            (x / c^6 + y / c + w / c^10, (z, y, z, z, z, z, x, z, z, z)),
            (x / c^6 - y / c, (z, -y, z, z, z, z, x)),
            (x / c^6 - y / c - w, (-w, -y, z, z, z, z, x)),
            (x / c^6 - y / c - w / c^10, (z, -y, z, z, z, z, x, z, z, z)),
            (x / c^6 - (y / c - w), (w, -y, z, z, z, z, x)),
            (x / c^6 - (y / c - w / c^10), (z, -y, z, z, z, z, x, z, z, z)),
            (-(x / c^6) - y / c, (z, -y, z, z, z, z, -x)),
            (-(x / c^6) - y / c - w, (-w, -y, z, z, z, z, -x)),
            (-(x / c^6) - y / c - w / c^10, (z, -y, z, z, z, z, -x, z, z, z)),
            ((x * c^2) / c^4 + y / c, (z, y, x)),
            ((x * c^2) / c^4 + y / c + w, (w, y, x)),
            ((x * c^2) / c^9 + y / c + w, (w, y, z, z, z, z, z, x)),
            (w + (x * c^2) / c^4 + y / c, (w, y, x)),
            (w + (x * c^2) / c^9 + y / c, (w, y, z, z, z, z, z, x)),
            (w + ((x * c^2) / c^4 + y / c), (w, y, x)),
            (w + ((x * c^2) / c^9 + y / c), (w, y, z, z, z, z, z, x)),
            (((x * c^2) / c^4 + y / c + w) / c^3, (z, z, z, w, y, x)),
            (((x * c^2) / c^4 + y / c + w) / c^5, (z, z, z, z, z, w, y, x)),
            (((x * c^2) / c^7 + y / c + w) / c^5, (z, z, z, z, z, w, y, z, z, z)),
            (v * (((x * c^2) / c^4 + y / c + w) / c^3), v .* (z, z, z, w, y, x)),
            ((((x * c^2) / c^4 + y / c + w) / c^3) * v, v .* (z, z, z, w, y, x)),
            ((((x * c^2) / c^4 + y / c + w) / c^3) / v, (z, z, z, w, y, x) .* (1 / v)),
        )
            @test expr.coeffs == expected
            @test pn_order(expr) == pn_order(pn)
            @test sum(expr) == sum(expected)
            @test eltype(expr) === T
        end

        # Can't make a PNExpression with positive exponents
        @test_throws ArgumentError x * c + y
        @test_throws ArgumentError x * c + y / c
        @test_throws ArgumentError x * c^2 + y
        @test_throws ArgumentError x * c^2 + y / c
        @test_throws ArgumentError y + x * c
        @test_throws ArgumentError y / c + x * c
        @test_throws ArgumentError y + x * c^2
        @test_throws ArgumentError y / c + x * c^2
        @test_throws ArgumentError x * c^2 + (y / c + z / c^2)
    end
end

@testitem "PNExpansion algebra" begin
    using Symbolics
    using PostNewtonian: PNExpansion

    for N1 ∈ 1:9
        for N2 ∈ 1:9
            for NMax ∈ max(N1, N2):(N1 + N2 + 3)
                @variables c⁻¹ x[1:N1] y[1:N2] z
                poly(e::PNExpansion) = sum(e[i] * c⁻¹^(i - 1) for i ∈ 1:length(e))
                eˣ = PNExpansion(tuple(x...), NMax)
                eʸ = PNExpansion(tuple(y...), NMax)

                # Test sums
                polysum = simplify(poly(eˣ + eʸ); expand=true)
                sumpoly = simplify(poly(eˣ) + poly(eʸ); expand=true)
                Δ = simplify(polysum - sumpoly; expand=true)
                @test iszero(Δ)
                @test_throws ArgumentError eˣ + PNExpansion(tuple(z, x...), NMax + 1)
                @test_throws ArgumentError PNExpansion(tuple(z, x...), NMax + 1) + eˣ

                # Test products
                polyprod = simplify(poly(eˣ * eʸ); expand=true)
                prodpoly = simplify(
                    substitute(
                        simplify(poly(eˣ) * poly(eʸ); expand=true),
                        Dict([c⁻¹^n => 0 for n ∈ NMax:(2NMax + 3)]),
                    );
                    expand=true,
                )
                Δ = simplify(polyprod - prodpoly; expand=true)
                @test iszero(Δ)
                @test_throws ArgumentError eˣ * PNExpansion(tuple(z, x...), NMax + 1)
                @test_throws ArgumentError PNExpansion(tuple(z, x...), NMax + 1) * eˣ
            end
        end
    end
end
