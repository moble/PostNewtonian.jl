@testset verbose=true "Aqua quality assurance tests" begin
    Aqua.test_all(PostNewtonian; ambiguities=false)
end
@testset verbose=true "ExplicitImports tests" begin
    @test ExplicitImports.check_no_implicit_imports(PostNewtonian) === nothing
    @test ExplicitImports.check_no_stale_explicit_imports(PostNewtonian) === nothing
    @test ExplicitImports.check_all_qualified_accesses_via_owners(PostNewtonian) === nothing
end
