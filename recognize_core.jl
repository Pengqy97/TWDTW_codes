#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# AUTHOR: Shen Ruoque
# VERSION: v2021.12.14
#=
使用面积卡阈值识别分类（核心函数）
=#

"""
    gridarea(Δφ, Δθ, θ)

caluculate grid's area (m²).

```math
Area = \\int_{\\varphi}^{\\varphi + \\Delta\\varphi} \\int_{\\theta}^{\\theta+\\Delta\\theta}
R \\cos\\theta \\mathrm{d}\\varphi R \\mathrm{d}\\theta
= R^2 \\Delta\\varphi(\\sin(\\theta + \\Delta\\theta) - \\sin\\theta)
```

Arguments:

- `Δφ`: resolution of longitude (degrees east)
- `Δθ`: resolution of latitude (degrees north)
- `θ`: latitude of the top of the grid (degree)
"""
function gridarea(Δφ, Δθ, θ) :: Float64
    R = 6371008.8
    abs(R ^ 2 * deg2rad(Δφ) * (sind(θ + Δθ) - sind(θ)))
end

"""
    recgnize(dist, ref, croparea, validmax=nothing; proj)

Recognize crop planting place. Returns the threshold and the classified map.

Arguments:

- dist: distance matrix, contains the similaraty of crop
- ref: georeference
- croparea: planting area of the crop (m²)
- validmax: the maximum posible distance value for the identified crop
- proj: projection (ALBERS or WGS84), string
"""
function recognize(dist, ref, croparea, validmax=nothing; proj, rev=false)
    if proj == "ALBERS"
        pixarea = abs(ref[:geotransform][2] * ref[:geotransform][6])
        pixnum = floor(Int, croparea / pixarea)
        sorted = sort(dist[:]; rev=rev)
        threshold = sorted[pixnum]
        classified = dist .≤ threshold
        threshold, classified
    elseif proj == "WGS84"
        _, Δφ, _, θ0, _, Δθ = ref[:geotransform]
        θ = [θ0 + (i - 1) * Δθ for i = 1:ref[:height]]
        pixareas = repeat(gridarea.(Δφ, Δθ, θ), 1, ref[:width])'
        minnum = floor.(Int, croparea / pixareas[1, end])
        maxnum = ceil.(Int, croparea / pixareas[1, 1])
        # 大于 validmax 的像元不参与排序，节省计算时间
        valid = isnothing(validmax) ? .!isnan.(dist) : dist .< validmax
        if isnothing(validmax)
            valid = .!isnan.(dist)
            if sum(valid) < maxnum error("cropland area not enough") end
        else
            valid = dist .< validmax
            while sum(valid) < maxnum
                validmax *= 2
                valid = dist .< validmax
            end
        end
        distvalid = dist[valid]
        pixareas = pixareas[valid]
        println("number of valid pixels: $(sum(valid))")
        println("maximum number of crop pixels: $maxnum")
        valid = nothing
        order = sortperm(distvalid; rev=rev)[1:maxnum]
        distpixarea = pixareas[order]
        pixareas = nothing
        println("sort")
        GC.gc()
        threshold = []
        accum = sum(distpixarea[1:minnum-1])
        for i = minnum:maxnum
            accum += distpixarea[i]
            if accum > croparea
                threshold = distvalid[order[i-1]]
                break
            end
        end
        threshold, dist .≤ threshold
    end
end
