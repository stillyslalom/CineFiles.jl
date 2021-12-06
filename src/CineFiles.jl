module CineFiles

using ColorTypes: N0f8, N4f12, Gray
using LRUCache
using StaticArrays

import Base: eltype, length, size, getindex, firstindex, lastindex, iterate
export CineFile

# using FileIO
# add_format(format"CINE", [0x43, 0x49], ".cine")

struct Time64
    fractions::UInt32
    seconds::UInt32
end


struct WBgain
	R::Float32
	B::Float32
end

struct IMfilter
	Dim::Int32
	Shifts::Int32
	Bias::Int32
	Coef::SMatrix{5,5,Int32,5*5}
end

struct Rect
	Left::Int32
	Top::Int32
	Right::Int32
	Bottom::Int32
end

struct TimeCode
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

struct SetupHeaderq
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
	DescriptionOld::SVector{121,Char}
	Mark::UInt16
	Length::UInt16
	Res13::UInt16
	SigOption::UInt16
	BinChannels::Int16
	SamplesPerImage::UInt8
	BinName::SMatrix{8,11,Char,8*11}
	AnaOption::UInt16
	AnaChannels::Int16
	Res6::UInt8
	AnaBoard::UInt8
	ChOption::SVector{8,Int16}
	AnaGain::SVector{8,Float32}
	AnaUnit::SMatrix{8,6,Char,8*6}
	AnaName::SMatrix{8,11,Char,8*11}
	lFirstImage::Int32
	dwImageCount::UInt32
	nQFactor::Int16
	wCineFileType::UInt16
	szCinePath::SMatrix{4,65,Char,4*65}
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
	bFlipH::SVector{4,Bool}
	bFlipV::SVector{4,Bool}
	Grid::UInt32
	FrameRate::UInt32
	Shutter::UInt32
	EDRShutter::UInt32
	PostTrigger::UInt32
	FrameDelay::UInt32
	bEnableColor::SVector{4,Bool}
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
	WBGain::SVector{4,WBgain}
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
	bStampTime::SVector{4,Bool}
	SoundDest::UInt32
	FRPSteps::UInt32
	FRPImgNr::SVector{16,Int32}
	FRPRate::SVector{16,UInt32}
	FRPExp::SVector{16,UInt32}
	MCCnt::Int32
	MCPercent::SVector{64,Float32}
	CICalib::UInt32
	CalibWidth::UInt32
	CalibHeight::UInt32
	CalibRate::UInt32
	CalibExp::UInt32
	CalibEDR::UInt32
	CalibTemp::UInt32
	HeadSerial::SVector{4,UInt32}
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
	Description::SVector{4096,Char}
	RisingEdge::SVector{4,Bool}
	FilterTime::UInt32
	LongReady::SVector{4,Bool}
	ShutterOff::SVector{4,Bool}
	Res4::SVector{16,UInt8}
	bMetaWB::SVector{4,Bool}
	Hue::Int32
	BlackLevel::Int32
	WhiteLevel::Int32
	LensDescription::SVector{256,Char}
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
	ToneLabel::SVector{256,Char}
	TonePoints::Int32
	fTone::SVector{64,Float32}
	UserMatrixLabel::SVector{256,Char}
	EnableMatrices::SVector{4,Bool}
	fUserMatrix::SVector{9,Float32}
	EnableCrop::SVector{4,Bool}
	CropRect::Rect
	EnableResample::SVector{4,Bool}
	ResampleWidth::UInt32
	ResampleHeight::UInt32
	fGain16_8::Float32
	FRPShape::SVector{16,UInt32}
	TrigTC::TimeCode
	fPbRate::Float32
	fTcRate::Float32
	CineName::SVector{256,Char}
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
	setup::SetupHeader
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

		seek(f, cine.OffSetup)
		setup = readstruct(f, SetupHeader)

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
    1 <= frameidx <= h.cine.ImageCount || BoundsError(h.dt, frameidx)
    seek(f, h.imglocs[frameidx])
    skip(f, read(f, UInt32) - 4)
    read!(f, h.tmp)

    # non-allocating version of rotl90
    ind1, ind2 = axes(h.tmp)
    n = first(ind2) + last(ind2)
    for i in axes(h.tmp, 1), j in ind2
        frame[n - j, i] = h.tmp[i, j]
    end
    return frame
end

readframe!(filename, frame, h, frameidx) = open(f -> readframe!(f, frame, h, frameidx), filename)
readframe(f, h, frameidx) = readframe!(f, similar(h.tmp, size(h.tmp, 2), size(h.tmp, 1)), h, frameidx)

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
