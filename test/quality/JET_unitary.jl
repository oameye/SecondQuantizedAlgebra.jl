using SecondQuantizedAlgebra
using JET
using Test
using Symbolics: @variables

@testset "JET Bogoliubov and affine substitution entry points" begin
    h = FockSpace(:fock)
    a = Destroy(h, :a)

    numeric = JET.@report_call target_modules = (SecondQuantizedAlgebra,) ignore_missing_comparison = true Bogoliubov(
        a, [1 0; 0 1],
    )
    @test isempty(JET.get_reports(numeric))

    complex_numeric = JET.@report_call target_modules = (SecondQuantizedAlgebra,) ignore_missing_comparison = true Bogoliubov(
        a, ComplexF64[im 0; 0 -im],
    )
    @test isempty(JET.get_reports(complex_numeric))

    @variables u::Number v::Number
    symbolic_matrix = [u v; conj(v) conj(u)]
    symbolic = JET.@report_call target_modules = (SecondQuantizedAlgebra,) ignore_missing_comparison = true Bogoliubov(
        a, symbolic_matrix,
    )
    @test isempty(JET.get_reports(symbolic))

    U = Bogoliubov(a, symbolic_matrix)
    resolved = JET.@report_call target_modules = (SecondQuantizedAlgebra,) ignore_missing_comparison = true substitute(
        U, Dict(u => 5 // 3, v => 4 // 3),
    )
    @test isempty(JET.get_reports(resolved))
end
