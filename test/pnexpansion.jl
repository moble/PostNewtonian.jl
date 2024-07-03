@testset verbose=true "PNExpansion" begin
    using PostNewtonian: PNExpansion

    for N1 ∈ 1:8
        for N2 ∈ 1:8
            for MaxN ∈ max(N1,N2):(N1+N2+3)
                @variables c⁻¹ x[1:N1] y[1:N2]
                poly(e::PNExpansion) = sum(e[i]*c⁻¹^(i-1) for i ∈ 1:length(e))
                eˣ = PNExpansion(tuple(x...), MaxN)
                eʸ = PNExpansion(tuple(y...), MaxN)

                # Test sums
                polysum = simplify(poly(eˣ + eʸ), expand=true)
                sumpoly = simplify(poly(eˣ)+poly(eʸ), expand=true)
                Δ = simplify(polysum - sumpoly, expand=true)
                @test iszero(Δ)
                @test_throws ErrorException eˣ + PNExpansion(tuple(x...), MaxN+1)
                @test_throws ErrorException PNExpansion(tuple(x...), MaxN+1) + eˣ

                # Test products
                polyprod = simplify(poly(eˣ * eʸ), expand=true)
                prodpoly = simplify(
                    substitute(
                        simplify(poly(eˣ)*poly(eʸ), expand=true),
                        Dict([c⁻¹^n=>0 for n ∈ MaxN:N1+N2+2])
                    ),
                    expand=true
                )
                Δ = simplify(polyprod - prodpoly, expand=true)
                @test iszero(Δ)
                @test_throws ErrorException eˣ * PNExpansion(tuple(x...), MaxN+1)
                @test_throws ErrorException PNExpansion(tuple(x...), MaxN+1) * eˣ

            end
        end
    end

    pn = rand(BBH; PNOrder=9//2)
    c = PNExpansionParameter(pn)

    # Test behavior of `c` as the basic PNTerm
    @test c.exp == -1
    @test c.coeff == 1
    @test c.coeff isa eltype(pn)
    @test (c^2).exp == -2
    @test (c^2).coeff == 1
    @test (c^2).coeff isa eltype(pn)
    @test (1.2 * c^2).exp == -2
    @test (1.2 * c^2).coeff == 1.2
    @test (1.2 * c^2).coeff isa eltype(pn)
    @test ((1.2 * c)^2).exp == -2
    @test ((1.2 * c)^2).coeff == 1.2^2
    @test ((1.2 * c)^2).coeff isa eltype(pn)
    @test (1.2 * c^2) / c^4 == 1.2 / c^2
    for (term, exp, coeff) ∈ (
        (c, -1, 1), (c^2, -2, 1), (1.2c^2, -2, 1.2), ((1.2c)^2, -2, 1.2^2),
        ((1.2c^2) / c^4, 2, 1.2), ((1.2c^2) / c^-4, -6, 1.2),
    )
        @test term.exp == exp
        @test term.coeff == coeff
        @test term.coeff isa eltype(pn)
    end

    # Test PNExpressions
    for (expr, expected) ∈ (
        1.2 / c^6 + 3.4/c, (0.0, 3.4, 0.0, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.0),
        1.2 / c^6 + 3.4/c + 5.6/c^10, (0.0, 3.4, 0.0, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.0),
        (1.2 * c^2) / c^4 + 3.4/c, (0.0, 3.4, 1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (1.2 * c^2) / c^4 + 3.4/c + 5.6, (5.6, 3.4, 1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (1.2 * c^2) / c^9 + 3.4/c + 5.6, (5.6, 3.4, 1.2, 0.0, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0),
        ((1.2 * c^2) / c^4 + 3.4/c + 5.6) / c^3, (0.0, 0.0, 0.0, 5.6, 3.4, 1.2, 0.0, 0.0, 0.0, 0.0),
        ((1.2 * c^2) / c^4 + 3.4/c + 5.6) / c^5, (0.0, 0.0, 0.0, 0.0, 0.0, 5.6, 3.4, 1.2, 0.0, 0.0),
        ((1.2 * c^2) / c^7 + 3.4/c + 5.6) / c^5, (0.0, 0.0, 0.0, 0.0, 0.0, 5.6, 3.4, 0.0, 0.0, 0.0),
    )
        @test expr.coeffs == expected
    end

    # Can't make a PNExpression with positive exponents
    @test_throws ErrorException 1.2 * c + 3.4
    @test_throws ErrorException 1.2 * c + 3.4 / c
    @test_throws ErrorException 1.2 * c^2 + 3.4
    @test_throws ErrorException 1.2 * c^2 + 3.4/c

end
