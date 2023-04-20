@testset verbose=true "truncated series monoid" begin
    using PostNewtonian: truncated_series_inverse, truncated_series_product,
        truncated_series_ratio
    @variables a[1:10] b[1:10] v

    @testset "$N" for N ∈ 0:6
        cₐ = collect(a[1:N+1])
        cᵦ = collect(b[1:N+1])
        A = sum(c*v^(i-1) for (i,c) ∈ enumerate(cₐ))
        A⁻¹ = sum(c*v^(i-1) for (i,c) ∈ enumerate(truncated_series_inverse(cₐ)))
        B = sum(c*v^(i-1) for (i,c) ∈ enumerate(cᵦ))
        truncate_orders(expr) = substitute(expr, Dict(v^n=>0 for n ∈ N+1:2length(a)))

        # Note that we have to engage in some chicanery with a[1] because Symbolics doesn't
        # deal well with fractions, and converts to Float64 too readily.
        full_product = simplify(a[1]^(N+1)*A*A⁻¹, expand=true)
        truncated_product = simplify(
            truncate_orders(full_product) / a[1]^(N+1),
            expand=true
        )
        @test iszero(truncated_product - 1)

        full_product = simplify(A*B, expand=true)
        truncated_product = simplify(truncate_orders(full_product), expand=true)
        full_product = truncated_series_product(cₐ, cᵦ, v)
        @test iszero(simplify(truncated_product - full_product, expand=true))

        @test iszero(simplify(truncated_series_ratio(cᵦ, cᵦ, v) - 1, expand=true))
    end

end
