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
                @test_throws eˣ + PNExpansion(tuple(x...), MaxN+1)
                @test_throws PNExpansion(tuple(x...), MaxN+1) + eˣ

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
                @test_throws eˣ * PNExpansion(tuple(x...), MaxN+1)
                @test_throws PNExpansion(tuple(x...), MaxN+1) * eˣ

            end
        end
    end
end
