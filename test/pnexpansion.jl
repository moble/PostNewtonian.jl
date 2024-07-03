@testset verbose=true "PNExpansion" begin
    using PostNewtonian: PNExpansion

    for N ∈ 1:9
        @variables c⁻¹ x[1:N] y[1:N] z
        poly(e::PNExpansion) = sum(e[i]*c⁻¹^(i-1) for i ∈ 1:length(e))
        eˣ = PNExpansion(tuple(x...))
        eʸ = PNExpansion(tuple(y...))

        # Test sums
        polysum = simplify(poly(eˣ + eʸ), expand=true)
        sumpoly = simplify(poly(eˣ)+poly(eʸ), expand=true)
        Δ = simplify(polysum - sumpoly, expand=true)
        @test iszero(Δ)
        @test_throws ErrorException eˣ + PNExpansion(tuple(z, x...))
        @test_throws ErrorException PNExpansion(tuple(z, x...)) + eˣ

        # Test products
        polyprod = simplify(poly(eˣ * eʸ), expand=true)
        prodpoly = simplify(
            substitute(
                simplify(poly(eˣ)*poly(eʸ), expand=true),
                Dict([c⁻¹^n=>0 for n ∈ N:2N+3])
            ),
            expand=true
        )
        Δ = simplify(polyprod - prodpoly, expand=true)
        @test iszero(Δ)
        @test_throws ErrorException eˣ * PNExpansion(tuple(z, x...))
        @test_throws ErrorException PNExpansion(tuple(z, x...)) * eˣ
    end

    for T ∈ [Float64, Float16, Double64]
        pn = rand(BBH; v=T(1)/5, PNOrder=9//2)
        c = PNExpansionParameter(pn)
        z = zero(T)
        x = T(6)/5
        y = T(17)/5
        w = T(28)/5
        v = T(27)/13

        # Test behavior of `c` as the basic PNTerm
        for (term, c⁻¹exp, coeff) ∈ (
            (c, -1, 1), (c^2, -2, 1), (x*c^2, -2, x), (c^2*x, -2, x), (c^2/x, -2, 1/x),
            ((x*c)^2, -2, x^2), ((x*c^2) / c^4, 2, x), ((x*c^2) / c^-4, -6, x),
            ((x/c^2) / c^4, 6, x), ((x/c^2) * (y/c^4), 6, x*y),
        )
            @test term.c⁻¹exp == c⁻¹exp
            @test term.coeff == coeff
            @test term.coeff isa eltype(pn)
            @test length(term) == 1
            @test eltype(term) === T
        end

        # Test PNExpressions
        @test_throws ArgumentError PNExpansion(())
        for (expr, expected) ∈ (
            (w - (x - y/c), (w-x, y, z, z, z, z, z, z, z, z)),
            (w - (x/c - y), (w+y, -x, z, z, z, z, z, z, z, z)),
            (x - y/c, (x, -y, z, z, z, z, z, z, z, z)),
            (x/c - y, (-y, x, z, z, z, z, z, z, z, z)),
            (x / c^6 + y/c, (z, y, z, z, z, z, x, z, z, z)),
            (x / c^6 + y/c + w/c^10, (z, y, z, z, z, z, x, z, z, z)),
            (x / c^6 - y/c, (z, -y, z, z, z, z, x, z, z, z)),
            (x / c^6 - y/c - w, (-w, -y, z, z, z, z, x, z, z, z)),
            (x / c^6 - y/c - w/c^10, (z, -y, z, z, z, z, x, z, z, z)),
            (x / c^6 - (y/c - w), (w, -y, z, z, z, z, x, z, z, z)),
            (x / c^6 - (y/c - w/c^10), (z, -y, z, z, z, z, x, z, z, z)),
            (-(x / c^6) - y/c, (z, -y, z, z, z, z, -x, z, z, z)),
            (-(x / c^6) - y/c - w, (-w, -y, z, z, z, z, -x, z, z, z)),
            (-(x / c^6) - y/c - w/c^10, (z, -y, z, z, z, z, -x, z, z, z)),
            ((x * c^2) / c^4 + y/c, (z, y, x, z, z, z, z, z, z, z)),
            ((x * c^2) / c^4 + y/c + w, (w, y, x, z, z, z, z, z, z, z)),
            ((x * c^2) / c^9 + y/c + w, (w, y, z, z, z, z, z, x, z, z)),
            (w + (x * c^2) / c^4 + y/c, (w, y, x, z, z, z, z, z, z, z)),
            (w + (x * c^2) / c^9 + y/c, (w, y, z, z, z, z, z, x, z, z)),
            (w + ((x * c^2) / c^4 + y/c), (w, y, x, z, z, z, z, z, z, z)),
            (w + ((x * c^2) / c^9 + y/c), (w, y, z, z, z, z, z, x, z, z)),
            (((x * c^2) / c^4 + y/c + w) / c^3, (z, z, z, w, y, x, z, z, z, z)),
            (((x * c^2) / c^4 + y/c + w) / c^5, (z, z, z, z, z, w, y, x, z, z)),
            (((x * c^2) / c^7 + y/c + w) / c^5, (z, z, z, z, z, w, y, z, z, z)),
            (v*(((x * c^2) / c^4 + y/c + w) / c^3), v.*(z, z, z, w, y, x, z, z, z, z)),
            ((((x * c^2) / c^4 + y/c + w) / c^3) * v, v.*(z, z, z, w, y, x, z, z, z, z)),
            ((((x * c^2) / c^4 + y/c + w) / c^3) / v, (z, z, z, w, y, x, z, z, z, z) .* (1/v)),
        )
            @test expr.coeffs == expected
            @test pn_order(expr) == pn_order(pn)
            @test sum(expr) == sum(expected)
            @test eltype(expr) === T
        end

        # Can't make a PNExpression with positive exponents
        @test_throws ArgumentError x*c + y
        @test_throws ArgumentError x*c + y / c
        @test_throws ArgumentError x*c^2 + y
        @test_throws ArgumentError x*c^2 + y/c
        @test_throws ArgumentError y + x*c
        @test_throws ArgumentError y / c + x*c
        @test_throws ArgumentError y + x*c^2
        @test_throws ArgumentError y / c + x*c^2
        @test_throws ArgumentError x*c^2 + (y / c + z / c^2)
    end
end
