module CineFiles

using ColorTypes: N0f8, N4f12, Gray
using LRUCache

import Base: eltype, length, size, getindex, firstindex, lastindex, iterate
export CineFile

# using FileIO
# add_format(format"CINE", [0x43, 0x49], ".cine")

struct Time64
    fractions::UInt32
    seconds::UInt32
end

struct CineFileHeader
    Type::UInt16
    HeaderSize::UInt16
    Compression::UInt16
    Version::UInt16
    FirstMovieImage::Int32
    TotalImageCount::UInt32
    FirstImageNum::Int32
    ImageCount::UInt32
    OffImageHeader::UInt32
    OffSetup::UInt32
    OffImageOffsets::UInt32
    TriggerTime::Time64
end

struct BitmapInfoHeader
    Size::UInt32
    Width::Int32
    Height::Int32
    Planes::UInt16
    BitCount::UInt16
    Compression::UInt32
    SizeImage::UInt32
    XPelsPerMeter::Int32
    YPelsPerMeter::Int32
    ClrUsed::UInt32
    ClrImportant::UInt32
    FPS::UInt16
end

readstruct(f, S::T) where {T} = only(read!(f, Vector{S}(undef, 1)))

struct CineHeader{T}
    cine::CineFileHeader
    bitmap::BitmapInfoHeader
    imglocs::Vector{Int64}
    imgoffset::Int
    dt::Vector{Float64}
    tmp::Array{T,2}
end

function CineHeader(fname)
    open(fname) do f
        read(f, UInt16) == UInt(18755) || error(basename(fname), " is not a .cine file")
        seek(f, 0)
        cine = readstruct(f, CineFileHeader)
        seek(f, cine.OffImageHeader)
        bitmap = readstruct(f, BitmapInfoHeader)

        bittype = bitmap.BitCount == 8 ? Gray{N0f8} : Gray{N4f12}
        cine.ImageCount > 0 || error("no images exist in file")

        seek(f, cine.OffImageOffsets)
        imglocs = read!(f, Array{Int64}(undef, cine.ImageCount))

        seekstart(f)
        imgoffset = 1
        while read(f, UInt16) != 1002
            imgoffset += 1
        end

        dt = zeros(cine.ImageCount)
        skip(f, 2)
        for i in eachindex(dt)
            fracstart = read(f, UInt32)
            secstart = read(f, UInt32)
            dt[i] = (secstart - cine.TriggerTime.seconds) + ((fracstart - cine.TriggerTime.fractions)/2^32)
        end

        tmp = Array{bittype, 2}(undef, bitmap.Width, bitmap.Height)
        return CineHeader{bittype}(cine, bitmap, imglocs, imgoffset, dt, tmp)
    end
end

Base.eltype(h::CineHeader{T}) where {T} = T

function readframe!(f::IO, frame, h, frameidx)
    frameidx <= h.cine.ImageCount || error("tried to access nonexistent frame $frameidx")
    seek(f, h.imglocs[frameidx])
    skip(f, read(f, UInt32) - 4)
    read!(f, h.tmp)
    frame .= rotl90(h.tmp)
end

readframe!(filename, frame, h, frameidx) = open(f -> readframe!(f, frame, h, frameidx), filename)
readframe(f, h, frameidx) = readframe!(f, rotl90(h.tmp), h, frameidx)

struct CineFile{T}
    path::String
    header::CineHeader{T}
    data::LRU{Int, Array{T,2}}
end

"""
    CineFile(filepath, cachelimit=0.25)

Load header information and create a frame cache for a Phantom .cine file.
    `cachelimit` sets the maximum size of the frame cache as a fraction of
    your system's free RAM. Indexing and iteration is supported for CineFiles.

"""
function CineFile(filepath, cachelimit=0.25)
    header = CineHeader(filepath)
    framesize = Base.summarysize(header.tmp)
    maxcachedframes = ceil(Int, cachelimit*Sys.free_memory() / framesize)
    data = LRU{Int, Array{eltype(header),2}}(maxsize = maxcachedframes)
    return CineFile(filepath, header, data)
end

Base.length(cf::CineFile) = length(cf.header.dt)
Base.eltype(cf::CineFile) = typeof(cf.header.tmp)
Base.size(cf::CineFile) = (length(cf), size(cf.header.tmp, 2), size(cf.header.tmp, 1))

function Base.getindex(cf::CineFile, idx::Int)
    get!(cf.data, idx) do
        readframe(cf.path, cf.header, idx)
    end
end
Base.firstindex(cf::CineFile) = 1
Base.lastindex(cf::CineFile) = length(cf)
Base.getindex(cf::CineFile, I) = [cf[i] for i in I]

Base.iterate(cf::CineFile, state=1) = state > length(cf) ? nothing : (cf[state], state + 1)

end
