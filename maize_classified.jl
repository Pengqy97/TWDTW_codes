#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# AUTHOR: Shen Ruoque
# VERSION: v2021.12.14
#=
识别玉米种植区
=#

using Distributed
addprocs(5)
@everywhere include("recognize_core.jl")
"""
    readgtiff(file)

A MATLAB-like way to read the GeoTiff file
return 3-D data array, and reference.
"""
function readgtiff(file)
    AG.read(file) do ds
        AG.read(ds)
    end, gtiffref(file)
end

fname = ARGS[1]
yrmin = parse(Int, ARGS[2])
path = ARGS[3]
province = ARGS[4]
city = ARGS[5]

maize_areas = []

yrmax = maize_areas[province] |> keys |> maximum

prefix0 = "distance-"
prefix = "classified-"

infile = joinpath(path, prefix0 * fname)
println(infile)
dist, ref = readgtiff(infile)
dist = dropdims(dist; dims=3)
maize_area = maize_areas[province][min(yrmin, yrmax)]
threshold, classified = recognize(dist, ref, 10000maize_area, 5; proj="WGS84", rev=false)
println(path)
outpath = joinpath(path, prefix * fname)
writegtiff(outpath, UInt8.(classified), ref)


rmprocs(workers())
