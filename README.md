# CineFiles

[![Build Status](https://github.com/stillyslalom/CineFiles.jl/workflows/CI/badge.svg)](https://github.com/stillyslalom/CineFiles.jl/actions)

CineFiles.jl reads `.cine` video files produced by [Phantom](https://www.phantomhighspeed.com/) high-speed cameras. It parses the binary header, extracts per-frame timestamps, and returns linearized grayscale image data with an LRU frame cache for efficient repeated access.

```julia
julia> using CineFiles

julia> cf = CineFile("test/data/8bpp.cine")
202-frame CineFile{Gray{N0f8}}(16, 128)

julia> cf[1]
16×128 Matrix{Gray{N0f8}}

julia> cf[end-2:end]
3-element Vector{Matrix{Gray{N0f8}}}
```

Supported pixel formats:
- 8-bit grayscale (uncompressed)
- 12-bit grayscale (uncompressed)
- 10-bit packed grayscale (Phantom lossy compression)

### Frame timestamps
Per-frame timestamps relative to the trigger event are stored in the header as a `Float64` vector, useful for synchronizing with external data or verifying frame rates.
```julia
julia> cf.header.dt[1:3]  # seconds relative to trigger
3-element Vector{Float64}:
 -0.009950000047683716
 -0.009900000095367432
 -0.009850000031292439

julia> 1 / (cf.header.dt[2] - cf.header.dt[1])  # frame rate
20000.000000000218
```

### Header metadata
Camera settings and acquisition parameters are accessible through the parsed `SetupHeader`, `CineFileHeader`, and `BitmapInfoHeader` structs, whose fields follow the naming conventions of the Phantom SDK C headers.
```julia
julia> cf.header.setup.FrameRate
20000

julia> cf.header.bitmap.Width
128

julia> cf.header.cine.ImageCount
202
```

### Frame cache
Frames are cached in an LRU cache sized as a fraction of free system RAM (default 25%). Set `cachelimit=0` to disable caching, or increase it for large batch processing.
```julia
julia> cf = CineFile("video.cine", 0.5)  # use up to 50% of free RAM for cache
```

### Iteration
`CineFile` supports `length`, `getindex`, `firstindex`, `lastindex`, and `iterate`, so standard iteration patterns work directly.
```julia
julia> length(cf)
202

julia> mean_intensity = sum(Float64.(gray.(frame)) for frame in cf) ./ length(cf);
```

### Installation
CineFiles.jl is registered in Julia's General package registry.
```
pkg> add CineFiles
```
