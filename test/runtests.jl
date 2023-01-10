using CineFiles
using ImageCore
using TiffImages
using Test
using Glob
using XMLDict: parse_xml
using FixedPointNumbers: nbitsfrac
using Aqua

# In order to test against arbitary cine files, export tif image of first frame with gamma = 1
# and xml file using PCC. Ensure the cine file, tif image and xml file have the same name and 
# place in either data or proprietary_data folders.

struct cine_test_files
    cine_path::String
    xml_path::String
    first_frame_path::String
end

function cine_test_files(cine_path::String)
    base_name = splitext(cine_path)[1]
    xml_path = base_name * ".xml"
    first_frame_path = base_name * ".tif"
    paths = (cine_path, xml_path, first_frame_path)
    exists = isfile.(paths)
    if all(exists)
        cine_test_files(paths...)
    else
        @error "Missing file(s) $(paths[[.~exists...]])"
    end
end

function compare_images(cine_file, tiff_image, frame_no=1)
    cine_frame = cine_file[frame_no]
    cine_frame_gamma = (cine_frame) .+ cine_file.header.setup.fOffset
    # Have not included gamma adjustment yet
    cine_frame_gamma[cine_frame_gamma.<0] .= 0
    return all(isapprox.(cine_frame_gamma, tiff_image; atol=0.01))
end

cine_file_paths = cine_test_files.(glob("*.cine", "data"))
append!(cine_file_paths, cine_test_files.(glob("*.cine", "proprietary_data")))

@testset "CineFiles.jl" begin
    @testset "Dependency internals" begin
        # nbitsfrac is unexported, so we'll test to make sure it doesn't break
        @test nbitsfrac(N0f8) == 8
        @test nbitsfrac(N4f12) == 12
        @test nbitsfrac(N0f16) == 16
    end

    @time cf8 = CineFile(joinpath("data", "8bpp.cine"))
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

    @test cf12[1] ≈ cf8[1]

    first_load_time = (@elapsed cf12[2])
    cached_load_time = (@elapsed cf12[2])
    @show first_load_time, cached_load_time
    @test first_load_time > cached_load_time

    @testset "Tiff compare first frame" begin
        for file in cine_file_paths
            @info file.cine_path
            cf = CineFile(file.cine_path)
            tiff_img = TiffImages.load(file.first_frame_path)
            @test compare_images(cf, tiff_img)
        end
    end

    @testset "XML compare camera version" begin
        for file in cine_file_paths
            @info file.cine_path
            cf = CineFile(file.cine_path)
            xml_data = parse_xml(read(file.xml_path, String))
            @test Int(cf.header.setup.CameraVersion) ==
                  parse(Int, xml_data["CameraSetup"]["CameraVersion"])
        end
    end
end

Aqua.test_all(CineFiles)
