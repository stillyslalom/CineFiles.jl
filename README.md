# CineFiles

Use `cf = CineFile("path/to/file.cine")` to create an object containing metadata and a frame cache for a grayscale Phantom .cine file. Individual frames can be loaded with, e.g. `cf[1]`, or ranges with `cf[1:3]` or `cf[end-5:end]`.
