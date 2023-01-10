module CineFiles

using ColorTypes: Color, N0f8, N6f10, N4f12, N0f16, gray, Gray, AbstractGray
using ColorVectorSpace
using LRUCache
using FixedPointNumbers: Normed, reinterpret, nbitsfrac

import Base: eltype, length, size, getindex, firstindex, lastindex, iterate, read, read!, setindex!
export CineFile

include("LUT.jl")

abstract type RawFrame{T} end
abstract type BinaryData end

struct AlArray{S,T,N} <: AbstractArray{T,N} where {S}
    data::Array{T,N}
    function AlArray{S,T}(data) where {S,T}
        if S isa Tuple
            size(data) == S || error("Array size mismatch")
        end
        if S isa Int
            length(data) == S || error("Array size mismatch")
        end
        new{S,T,length(S)}(data)
    end
end

size(A::AlArray) = size(A.data)
size(::Type{AlArray{S,T}}) where {S,T} = Tuple(S)
getindex(A::AlArray{S,T,N}, inds::Vararg{Int,N}) where {S,T,N} = A.data[inds...]
setindex!(A::AlArray{S,T,N}, val, inds::Vararg{Int,N}) where {S,T,N} = A.data[inds...] = val

read(f::IOStream, A::Type{AlArray{S,T}}) where {S,T} = A(read!(f, Array{T,length(S)}(undef, S)))

rawcounts(v::AbstractGray) = reinterpret(gray(v))

struct Levels{T}
    white::T
    black::T
end

function Levels{T}(W::S, B::S) where {T<:AbstractGray,S<:Int32}
    scale = rawcounts(T(1))
    Levels{T}(T(W / scale), T(B / scale))
end

struct Unpacked{T} <: RawFrame{T}
    width::Int
    height::Int
    tmp::Array{T,2}
    unpacked::Array{T,2}
    levels::Levels{T}
    Unpacked{T}(width, height, whiteL, blackL) where {T<:Color} = new{T}(
        width,
        height,
        Array{T,2}(undef, width, height),
        Array{T,2}(undef, width, height),
        Levels{T}(whiteL, blackL),
    )
end

struct Packed{T,S} <: RawFrame{T}
    width::Int
    height::Int
    tmp::Array{S,2}
    unpacked::Array{T,2}
    pack::Vector{UInt8}
    levels::Levels{T}
end

function Packed{T,S}(width, height, whiteL, blackL) where {T<:Gray,S<:Gray}
    true_bits = nbitsfrac(eltype(T)) * width * height
    true_bytes = ceil(Int, true_bits / 8)
    Packed{T,S}(
        width,
        height,
        Array{S,2}(undef, width, height),
        Array{T,2}(undef, width, height),
        Vector{UInt8}(undef, true_bytes),
        Levels{T}(whiteL, blackL),
    )
end

function Packed{T}(width, height, whiteL, blackL) where {T<:Gray}
    true_bits = nbitsfrac(eltype(T)) * width * height
    true_bytes = ceil(Int, true_bits / 8)
    Packed{T,T}(
        width,
        height,
        Array{T,2}(undef, width, height),
        Array{T,2}(undef, width, height),
        Vector{UInt8}(undef, true_bytes),
        Levels{T}(whiteL, blackL),
    )
end

struct Time64 <: BinaryData
    fractions::UInt32
    seconds::UInt32
end

struct WBgain <: BinaryData
    R::Float32
    B::Float32
end

struct IMfilter <: BinaryData
    Dim::Int32
    Shifts::Int32
    Bias::Int32
    Coef::AlArray{(5, 5),Int32}
end

struct Rect <: BinaryData
    Left::Int32
    Top::Int32
    Right::Int32
    Bottom::Int32
end

struct TimeCode <: BinaryData
    framesU::UInt8
    framesT::UInt8
    dropFrameFlag::UInt8
    colorFrameFlag::UInt8
    secondsU::UInt8
    secondsT::UInt8
    flag1::UInt8
    minutesU::UInt8
    minutesT::UInt8
    flag2::UInt8
    hoursU::UInt8
    hoursT::UInt8
    flag3::UInt8
    flag4::UInt8
    userBitData::UInt32
end

struct SetupHeader <: BinaryData
    FrameRate16::UInt16
    Shutter16::UInt16
    PostTrigger16::UInt16
    FrameDelay16::UInt16
    AspectRatio::UInt16
    Res7::UInt16
    Res8::UInt16
    Res9::UInt8
    Res10::UInt8
    Res11::UInt8
    TrigFrame::UInt8
    Res12::UInt8
    DescriptionOld::AlArray{121,Int8}
    Mark::UInt16
    Length::UInt16
    Res13::UInt16
    SigOption::UInt16
    BinChannels::Int16
    SamplesPerImage::UInt8
    BinName::AlArray{(8, 11),Int8}
    AnaOption::UInt16
    AnaChannels::Int16
    Res6::UInt8
    AnaBoard::UInt8
    ChOption::AlArray{8,Int16}
    AnaGain::AlArray{8,Float32}
    AnaUnit::AlArray{(8, 6),Int8}
    AnaName::AlArray{(8, 11),Int8}
    lFirstImage::Int32
    dwImageCount::UInt32
    nQFactor::Int16
    wCineFileType::UInt16
    szCinePath::AlArray{(4, 65),Int8}
    Res14::UInt16
    Res15::UInt8
    Res16::UInt8
    Res17::UInt16
    Res18::Float64
    Res19::Float64
    Res20::UInt16
    Res1::Int32
    Res2::Int32
    Res3::Int32
    ImWidth::UInt16
    ImHeight::UInt16
    EDRShutter16::UInt16
    Serial::UInt32
    Saturation::Int32
    Res5::UInt8
    AutoExposure::UInt32
    bFlipH::Int32
    bFlipV::Int32
    Grid::UInt32
    FrameRate::UInt32
    Shutter::UInt32
    EDRShutter::UInt32
    PostTrigger::UInt32
    FrameDelay::UInt32
    bEnableColor::Int32
    CameraVersion::UInt32
    FirmwareVersion::UInt32
    SoftwareVersion::UInt32
    RecordingTimeZone::Int32
    CFA::UInt32
    Bright::Int32
    Contrast::Int32
    Gamma::Int32
    Res21::UInt32
    AutoExpLevel::UInt32
    AutoExpSpeed::UInt32
    AutoExpRect::Rect
    WBGain::AlArray{4,WBgain}
    Rotate::Int32
    WBView::WBgain
    RealBPP::UInt32
    Conv8Min::UInt32
    Conv8Max::UInt32
    FilterCode::Int32
    FilterParam::Int32
    UF::IMfilter
    BlackCalSVer::UInt32
    WhiteCalSVer::UInt32
    GrayCalSVer::UInt32
    bStampTime::Int32
    SoundDest::UInt32
    FRPSteps::UInt32
    FRPImgNr::AlArray{16,Int32}
    FRPRate::AlArray{16,UInt32}
    FRPExp::AlArray{16,UInt32}
    MCCnt::Int32
    MCPercent::AlArray{64,Float32}
    CICalib::UInt32
    CalibWidth::UInt32
    CalibHeight::UInt32
    CalibRate::UInt32
    CalibExp::UInt32
    CalibEDR::UInt32
    CalibTemp::UInt32
    HeadSerial::AlArray{4,UInt32}
    RangeCode::UInt32
    RangeSize::UInt32
    Decimation::UInt32
    MasterSerial::UInt32
    Sensor::UInt32
    ShutterNs::UInt32
    EDRShutterNs::UInt32
    FrameDelayNs::UInt32
    ImPosXAcq::UInt32
    ImPosYAcq::UInt32
    ImWidthAcq::UInt32
    ImHeightAcq::UInt32
    Description::AlArray{4096,Int8}
    RisingEdge::Int32
    FilterTime::UInt32
    LongReady::Int32
    ShutterOff::Int32
    Res4::AlArray{16,UInt8}
    bMetaWB::Int32
    Hue::Int32
    BlackLevel::Int32
    WhiteLevel::Int32
    LensDescription::AlArray{256,Int8}
    LensAperture::Float32
    LensFocusDistance::Float32
    LensFocalLength::Float32
    fOffset::Float32
    fGain::Float32
    fSaturation::Float32
    fHue::Float32
    fGamma::Float32
    fGammaR::Float32
    fGammaB::Float32
    fFlare::Float32
    fPedestalR::Float32
    fPedestalG::Float32
    fPedestalB::Float32
    fChroma::Float32
    ToneLabel::AlArray{256,Int8}
    TonePoints::Int32
    fTone::AlArray{64,Float32}
    UserMatrixLabel::AlArray{256,Int8}
    EnableMatrices::Int32
    fUserMatrix::AlArray{9,Float32}
    EnableCrop::Int32
    CropRect::Rect
    EnableResample::Int32
    ResampleWidth::UInt32
    ResampleHeight::UInt32
    fGain16_8::Float32
    FRPShape::AlArray{16,UInt32}
    TrigTC::TimeCode
    fPbRate::Float32
    fTcRate::Float32
    CineName::AlArray{256,UInt8}
end

struct CineFileHeader <: BinaryData
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

struct BitmapInfoHeader <: BinaryData
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

struct CineHeader{T, R}
    cine::CineFileHeader
    bitmap::BitmapInfoHeader
    setup::SetupHeader
    imglocs::Vector{Int64}
    imgoffset::Int
    dt::Vector{Float64}
    raw::R
end

struct CineFile{T,R<:RawFrame}
    path::String
    header::CineHeader{T,R}
    data::LRU{Int,Array{T,2}}
end

function unpack!(pack::Vector{UInt8}, unpacked::Array{Gray{N6f10}})
    pack16 = UInt16.(pack)
    indx = CartesianIndices(unpacked)[:, end:-1:1]
    unpacked[indx[1:4:end]] =
        reinterpret.(N6f10, (pack16[1:5:end] .<< 2) .| (pack16[2:5:end] .>> 6))
    unpacked[indx[2:4:end]] =
        reinterpret.(
            N6f10,
            ((pack16[2:5:end] .& 0b00111111) .<< 4) .| (pack16[3:5:end] .>> 4),
        )
    unpacked[indx[3:4:end]] =
        reinterpret.(
            N6f10,
            ((pack16[3:5:end] .& 0b00001111) .<< 6) .| (pack16[4:5:end] .>> 2),
        )
    unpacked[indx[4:4:end]] =
        reinterpret.(N6f10, ((pack16[4:5:end] .& 0b00000011) .<< 8) .| pack16[5:5:end])
end

function unpack!(pack::Vector{UInt8}, unpacked::Array{Gray{N4f12}})
    @warn "unpacking 12 bit images untested"
    pack16 = UInt16.(pack)
    indx = CartesianIndices(unpacked)[:, end:-1:1]
    unpacked[indx[1:2:end]] =
        reinterpret.(N6f10, (pack16[1:3:end] .<< 4) .| (pack16[2:3:end] .>> 4))
    unpacked[indx[2:2:end]] =
        reinterpret.(N6f10, ((pack16[2:3:end] .& 0b00001111) .<< 8) .| (pack16[3:3:end]))
end

function linearize!(raw_data::R) where {R<:RawFrame}
    Wlevel = raw_data.levels.white
    Blevel = raw_data.levels.black
    raw_data.tmp .= raw_data.unpacked
    clamp!(raw_data.tmp, Blevel, Wlevel)
    @inbounds for idx in eachindex(raw_data.tmp)
        raw_data.tmp[idx] = (raw_data.tmp[idx] - Blevel) / (Wlevel - Blevel)
    end
end

function linearize!(raw_data::Packed{Gray{N6f10}})
    Blevel, Wlevel = Gray{N4f12}.((64, 4064) ./ (2^12))
    raw_data.tmp .= lookup.(raw_data.unpacked, Ref(CINE_LUT))
    raw_data.tmp[raw_data.tmp .> Wlevel] .= Wlevel
    raw_data.tmp[raw_data.tmp .< Blevel] .= Blevel
    raw_data.tmp .= (raw_data.tmp .- Blevel) ./ (Wlevel - Blevel)
end

function CineHeader(fname)
    open(fname) do f
        read(f, UInt16) == UInt(18755) || error(basename(fname), " is not a .cine file")
        seek(f, 0)
        cine = read(f, CineFileHeader)
        seek(f, cine.OffImageHeader)
        bitmap = read(f, BitmapInfoHeader)

        seek(f, cine.OffSetup)
        setup = read(f, SetupHeader)

        bitsized = (bitmap.BitCount == 8) || (bitmap.BitCount == 24) ? UInt8 : UInt16

        if (bitmap.BitCount == 8) || (bitmap.BitCount == 16)
            bittype = Gray{Normed{bitsized,Int(setup.RealBPP)}}
        elseif (bitmap.BitCount == 24) || (bitmap.BitCount == 36)
            bittype = BGR{Normed{bitsized,Int(setup.RealBPP)}}
        end

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
            dt[i] =
                (secstart - cine.TriggerTime.seconds) +
                ((fracstart - cine.TriggerTime.fractions) / 2^32)
        end

        # Packing options
        if bitmap.Compression == 0
            raw = Unpacked{bittype}(
                bitmap.Width,
                bitmap.Height,
                setup.WhiteLevel,
                setup.BlackLevel,
            )
            pixeltype = bittype
        elseif bitmap.Compression == 256
            raw = Packed{bittype,Gray{N6f10}}(
                bitmap.Width,
                bitmap.Height,
                setup.WhiteLevel,
                setup.BlackLevel,
            )
            pixeltype = Gray{N4f12}
        elseif bitmap.Compression == 1024
            raw = Packed{bittype,Gray{N4f12}}(
                bitmap.Width,
                bitmap.Height,
                setup.WhiteLevel,
                setup.BlackLevel,
            )
            pixeltype = bittype
        end
        return CineHeader{pixeltype,typeof(raw)}(cine, bitmap, setup, imglocs, imgoffset, dt, raw)
    end
end

function readframe!(f::IO, frame, h, frameidx)
    1 <= frameidx <= h.cine.ImageCount || BoundsError(h.dt, frameidx)
    seek(f, h.imglocs[frameidx])
    skip(f, read(f, UInt32) - 4)
    read!(f, h.raw)
    linearize!(h.raw)
    # non-allocating version of rotl90
    ind1, ind2 = axes(h.raw.tmp)
    n = first(ind2) + last(ind2)
    for i in axes(h.raw.tmp, 1), j in ind2
        frame[n-j, i] = h.raw.tmp[i, j]
    end
    return frame
end

function read!(f::IOStream, frame::Packed)
    read!(f, frame.pack)
    unpack!(frame.pack, frame.unpacked)
end

read!(f::IOStream, frame::Unpacked) = read!(f, frame.unpacked)
# TODO: implement with lower compilation overhead
read(f::IOStream, S::Type{T}) where {T<:BinaryData} = S((read(f, FT) for FT in fieldtypes(T))...)# in field.(Ref(f), S.types)...)

readframe!(filename, frame, h, frameidx) =
    open(f -> readframe!(f, frame, h, frameidx), filename)
readframe(f, h, frameidx) =
    readframe!(f, similar(h.raw.tmp, size(h.raw.tmp, 2), size(h.raw.tmp, 1)), h, frameidx)

"""
    CineFile(filepath, cachelimit=0.25)

Load header information and create a frame cache for a Phantom .cine file.
    `cachelimit` sets the maximum size of the frame cache as a fraction of
    your system's free RAM. Indexing and iteration is supported for CineFiles.

"""
function CineFile(filepath, cachelimit=0.25)
    header = CineHeader(filepath)
    framesize = Base.summarysize(header.raw.tmp)
    maxcachedframes = ceil(Int, cachelimit * Sys.free_memory() / framesize)
    data = LRU{Int,Array{eltype(header),2}}(maxsize=maxcachedframes)
    return CineFile(filepath, header, data)
end

Base.eltype(h::CineHeader{T}) where {T} = T
Base.length(cf::CineFile) = length(cf.header.dt)
Base.eltype(cf::CineFile) = typeof(cf.header.raw.tmp)
Base.size(cf::CineFile) = (length(cf), size(cf.header.raw.tmp, 2), size(cf.header.raw.tmp, 1))

function Base.getindex(cf::CineFile, idx::Int)
    get!(cf.data, idx) do
        readframe(cf.path, cf.header, idx)
    end
end

Base.firstindex(cf::CineFile) = 1
Base.lastindex(cf::CineFile) = length(cf)
Base.getindex(cf::CineFile, I) = [cf[i] for i in I]
Base.iterate(cf::CineFile, state=1) = state > length(cf) ? nothing : (cf[state], state + 1)

# precompile as the final step of the module definition:
if ccall(:jl_generating_output, Cint, ()) == 1   # if we're precompiling the package
    let
        datadir = joinpath(splitpath(@__DIR__)[1:end-1]..., "test", "data")
        s = Gray{Float32}(0.0)
        for fmt in ("8bpp", "12bpp", "packed_10")
            collect(CineFile(joinpath(datadir, "$fmt.cine"), 0))
        end
    end
end

end
