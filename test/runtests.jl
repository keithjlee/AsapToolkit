using Test
using Asap
using AsapToolkit

#=
AsapToolkit v1.0 test suite.

The generators, Geo extraction, FDM translations, and model utilities moved
into Asap itself (with their tests — see Asap's test/newcore/test_generation.jl).
What remains here: the AISC steel-section database, JSON/Grasshopper IO,
and the AsapSections polygon-geometry subpackage.
=#

#a basic frame section (kN-m-ish units)
const mat = Material(200e6, 77e6, 80.0, 0.3)
const frame_section = Section(mat, 1e-2, 1e-4, 5e-5, 1e-6)
const truss_section = Section(mat, 1e-3)

@testset "AsapToolkit" begin

    @testset "SteelSections" begin
        w = W(String(first(Wnames)))
        @test w.A > 0

        fsec = toASAPframe(w, 200.0, 77.0)
        @test fsec isa Section
        @test EA(fsec) > 0
        @test EIx(fsec) > 0
        @test GJ(fsec) > 0

        tsec = toASAPtruss(w, 200.0)
        @test tsec isa Section
        @test EA(tsec) > 0
    end

    @testset "IO" begin
        #generators now come from Asap
        truss_gen = Pratt2D(12.0, 6, 1.5, truss_section)

        topo = AsapToolkit.topologize(truss_gen.model)
        @test topo isa Dict
        @test length(topo["iStart"]) == length(truss_gen.model.elements)

        dir = mktempdir()
        fname = joinpath(dir, "ghtest")
        GHsave(truss_gen.model, fname)
        @test isfile(fname * ".json")

        frame_gen = GridFrame(6.0, 4, 6.0, 4, frame_section)
        fname2 = joinpath(dir, "ghframe")
        GHsave(frame_gen.model, fname2)
        @test isfile(fname2 * ".json")
    end

    @testset "AsapSections" begin
        #unit square: area, centroid
        pts = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
        sec = SolidSection(hcat(pts...))
        @test sec.area ≈ 1.0
        @test sec.centroid ≈ [0.5, 0.5]

        props = SectionProperties(sec)
        @test props.area ≈ 1.0
    end

end
