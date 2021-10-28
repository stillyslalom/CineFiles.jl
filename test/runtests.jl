using CineFiles
using ImageCore
using TiffImages
using Test

@testset "CineFiles.jl" begin
    cf8 = CineFile(joinpath("data", "8bpp.cine"))
    @testset "8 bit grayscale" begin
        @test length(cf8) == 202
        @test size(cf8[1]) == (16, 128)
        @test eltype(cf8) == Matrix{Gray{N0f8}}
    end

    cf12 = CineFile(joinpath("data", "12bpp.cine"))
    @testset "12 bit grayscale" begin
        @test length(cf12) == 202
        @test size(cf12[1]) == (16, 128)
        @test eltype(cf12) == Matrix{Gray{N4f12}}
    end

    @test Float32.(cf12[1]) ≈ Float32.(cf8[1]) atol=0.1

    first_load_time = (@elapsed cf12[2])
    cached_load_time = (@elapsed cf12[2])
    @show first_load_time, cached_load_time
    @test first_load_time > cached_load_time

    # TODO: add comparison to tiff. 
    # TiffImages.jl cannot currently read the Phantom-exported tiff format.
end
