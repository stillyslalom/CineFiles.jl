# CineFiles

Use `cf = CineFile("path/to/file.cine")` to create an object containing metadata and a frame cache for a video file in .cine format produced by a [Phantom](https://www.phantomhighspeed.com/) high-speed camera. Individual frames can be loaded with, e.g. `cf[1]`, or ranges with `cf[1:3]` or `cf[end-5:end]`.

## Installation
`CineFiles` is included in Julia's `General` package registry and can be added to the current environment with `Pkg.add("CineFiles")`.
