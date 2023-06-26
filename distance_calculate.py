#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   distance_calculate.py
@Create Time    :   2022/02/25 14:04:04
@Author  :   pengqy 
'''
# Purpose: 按玉米生长季识别种植位置，计算每年的距离值。
# %%
from os.path import join
from genericpath import isfile
import rasterio, glob
from osgeo import gdal
import numpy as np
import os, time, ctypes
import rasterio.mask
import multiprocessing as mp
from datetime import timedelta
from twdtw import tw_dtw, dtw_core

#%%
julia = "/usr/local/julia/julia-1.6.0/bin/julia"
def readgtiff(file):
    with rasterio.open(file) as ds:
        data = ds.read()
        ref = {
            "width": ds.width,
            "height": ds.height,
            "transform": ds.transform,
            "crs": ds.crs
        }
    return data, ref

def writegtiff(file, data, ref):
    with rasterio.open(
        file, "w", driver="GTiff",
        width=ref["width"], height=ref["height"], count=1,
        dtype=data.dtype,
        crs=ref["crs"], transform=ref["transform"],
        compress="DEFLATE"
    ) as ds:
        ds.write(data[:, :], 1)

def f_summer_maize(pix1, threshold=500, n=9, f=0.25):
    """夏玉米权重"""
    return 1 + f * np.max([0, np.sum(pix1 > threshold) - n])

def get_dist(i):
    scale = np.float32(1000)
    data = np.frombuffer(data_shared, dtype=ctypes.c_int16)
    data = data.reshape(len(t), data.size // (len(t)))
    pix = data[:, i]
    tw0 = tws
    dist = (
        tw_dtw.single16(pix, std, tw0, scale) 
        * f_summer_maize(pix)
    )
    return dist
#%%
province = "Tianjin"
print(province)    
outpath = ""
if not os.path.exists(outpath):
    os.mkdir(outpath)
inputpath = ""
version = "1"
d_time = {
    "Tianjin" : [169, 297], 
}
d_start = d_time[province][0]; d_end = d_time[province][1]
d_start_n = int((d_start - 1) / 8); d_end_n = int((d_end - 1) / 8) + 1
std0 = {
    "Tianjin" : np.array([
        0.185142857, 0.175238095, 0.166904762, 0.158690476, 0.151904762, 
        0.148380952, 0.150071429, 0.157452381, 0.173952381, 0.203928571, 
        0.250904762, 0.314190476, 0.390809524, 0.474309524, 0.55452381, 
        0.618380952, 0.650904762, 0.65997619, 0.653333333, 0.641452381, 
        0.633952381, 0.63697619, 0.652357143, 0.679214286, 0.713166667, 
        0.74897619, 0.782190476, 0.80852381, 0.825690476, 0.830071429, 
        0.819380952, 0.791547619, 0.746071429, 0.682285714, 0.600928571, 
        0.5125, 0.428119048, 0.356214286, 0.303285714, 0.268547619, 
        0.246047619, 0.234452381, 0.23197619, 0.235785714, 0.243285714, 
        0.249571429
    ], dtype=np.float32), ##NDVI-Tianjin
}
std = std0[province][d_start_n:d_end_n]
for yr in range(2001,2020+1):

    print("caculate distance")

    filname = "distance-"
    outfile = os.path.join(outpath, filname)

    inpath_flist = glob.glob(inputpath + f"/*Y{yr}*{province}.tif")
    Province_CROP_ndvi = np.array([gdal.Open(tif_file).ReadAsArray() for tif_file in inpath_flist])
    _, ref0 = readgtiff(inpath_flist[0])

    valid = np.any((Province_CROP_ndvi > 300), axis=0)
    index = np.nonzero(valid)
    validnum = len(index[0])
    dists = ~valid * np.float32(np.nan)

    Province_CROP_ndvi = Province_CROP_ndvi[:, index[0], index[1]]
    data_shared = Province_CROP_ndvi.copy().flatten()

    '''时间矩阵'''
    tws = {}
    t = np.arange(d_start, d_end+1, 8, dtype=np.int64)
    t_std = np.arange(d_start, d_end+1, 8, dtype=np.int64)
    tws = dtw_core.tw_single64(t, t_std, 0.1, 50.0)

    '''并行计算距离值'''
    if not os.path.isfile(outfile):
        pools = mp.Pool(processes=80)
        start_time = time.time()
        print("start: ", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

        distance_list = pools.map(get_dist, range(validnum))

        end_time = time.time()
        print("end: ", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        print("total time: ", timedelta(seconds=end_time-start_time))

        dists[index[0], index[1]] = distance_list
        writegtiff(outfile, dists, ref0)
        pools.close()
        pools.join()

    print("classifing")


