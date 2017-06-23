from __future__ import division
from __future__ import print_function
'''
Zonal Statistics
Vector-Raster Analysis

Copyright 2013 Matthew Perry

All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of the software nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Usage:
  zonal_stats.py VECTOR RASTER
  zonal_stats.py -h | --help
  zonal_stats.py --version

Options:
  -h --help     Show this screen.
  --version     Show version.
'''
import os,sys
from osgeo import gdal, ogr

from osgeo.gdalconst import *
import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import zoom
import csv
from math import degrees, atan
import matplotlib.pyplot as plt
gdal.PushErrorHandler('CPLQuietErrorHandler')
 
 
def bbox_to_pixel_offsets(gt, bbox, rsize):
    originX = gt[0]
    originY = gt[3]
    pixel_width = gt[1]
    pixel_height = gt[5]
    x1 = int((bbox[0] - originX) / pixel_width)
    x2 = int((bbox[1] - originX) / pixel_width) + 1
 
    y1 = int((bbox[3] - originY) / pixel_height)
    y2 = int((bbox[2] - originY) / pixel_height) + 1
    
    # "Clip" the geometry bounds to the overall raster bounding box
    # This should avoid any rasterIO errors for partially overlapping polys
    if x1 < 0:
        x1 = 0
    if x2 > rsize[0]:
        x2 = rsize[0]
    if y1 < 0:
        y1 = 0
    if y2 > rsize[1]:
        y2 = rsize[1]

    xsize = x2 - x1
    ysize = y2 - y1

    return (x1, y1, xsize, ysize)
 
 
def zonal_stats(gdb, fc, fldname, ras, buffDist, filenm, nodata_value=None):
    #Open raster
    rds = gdal.Open(ras, GA_ReadOnly)
    assert(rds)
    rb = rds.GetRasterBand(1)
    rgt = rds.GetGeoTransform()
    # print(rgt)
    rsize = (rds.RasterXSize, rds.RasterYSize)
   #Get raster cell size
   #this works in case cell size is not square
   #for square rasters len1 = len2
    dx1 = rgt[1]
    dy1 = rgt[2]
    dx2 = rgt[4]
    dy2 = rgt[5]
    len1 = np.sqrt(dx1**2 + dy1 **2)
    len2 = np.sqrt(dx2**2 + dy2 **2)
    ras_area = len1 * len2
    #if needed, angle of rotation (rotated rasters)
    theta = np.arctan(dy1/dx1)

   #Establish NoData value -- note: can be provided as input
    if nodata_value:
        nodata_value = float(nodata_value)
        rb.SetNoDataValue(nodata_value)
    else:
        nodata_value = rb.GetNoDataValue()
    print('Raster NODATA Value: ', nodata_value)
    #Open feature class
    vds = ogr.Open(gdb, GA_ReadOnly)  
    vlyr = vds.GetLayerByName(fc)


    #Create memory drivers to hold arrays
    mem_drv = ogr.GetDriverByName('Memory')
    driver = gdal.GetDriverByName('MEM')
 
    # Loop through vectors
    statDict = {}
    lenid = len(vlyr)
    for i, feat in enumerate(vlyr):
        fldid = feat.GetField(fldname)
        print('{} of {}, staid: {}'.format(i+1, lenid, fldid))
        #Buffer well points, using buffDist input (in meters)
        geom = feat.GetGeometryRef()
        buff = geom.Buffer(buffDist) 
        vec_area = buff.GetArea()
        # print('Ratio Buffer to Raster: ', vec_area/ras_area)
        src_offset = bbox_to_pixel_offsets(rgt, buff.GetEnvelope(), rsize)
        if src_offset[2] <= 0 or src_offset[3] <=0:
            #if point falls outside raster grid, include nodata as zone analysis
            statDict[feat.GetField(fldname)] = nodata_value
        else:    
            if vec_area/ras_area >= 4:
                zooms = 1
            else:
                zooms = int(np.sqrt(ras_area/(vec_area/4)))
            
            src_array = rb.ReadAsArray(*src_offset)
            # print(src_array)
            # plt.imshow(src_array)
            # plt.show()

            #Calculate new geotransform of the feature subset
            new_gt = ((rgt[0] + (src_offset[0] * rgt[1])),
                rgt[1]/zooms,
                0.0,
                (rgt[3] + (src_offset[1] * rgt[5])),
                0.0,
                rgt[5]/zooms)

            #Create a temporary vector layer in memory
            mem_ds = mem_drv.CreateDataSource('out')
            mem_layer = mem_ds.CreateLayer('poly', None, ogr.wkbPolygon)
            mem_poly = ogr.Feature(mem_layer.GetLayerDefn())
            mem_poly.SetGeometryDirectly(buff)
            mem_layer.CreateFeature(mem_poly)
            
            #Rasterize the polygon
            rvds = driver.Create('', src_offset[2]*zooms, src_offset[3]*zooms, 1, gdal.GDT_Byte)
            rvds.SetGeoTransform(new_gt)
            gdal.RasterizeLayer(rvds, [1], mem_layer, burn_values=[1])
            rv_array = rvds.ReadAsArray()
            # print(rv_array)
            # plt.imshow(rv_array)
            # plt.show()

            #Resample the raster (only changes when zooms not 1)
            src_re = zoom(src_array, zooms, order = 0)
            # print(src_re)
            # plt.imshow(src_re)
            # plt.show()
            
            # Mask the source data array with our current feature
            # we take the logical_not to flip 0<->1 to get the correct mask effect
            # we also mask out nodata values explictly
            masked = np.ma.MaskedArray(src_re,
                mask=np.logical_or(
                    src_re == nodata_value,
                    np.logical_not(rv_array)))
            # print(masked)
            # plt.imshow(masked)
            # plt.show()
            #NOTE: Have tried to add mode functionality, but needs work
            # print(stats.mode(masked))

            feature_stats = {
                'min': float(masked.min()),
                'mean': float(masked.mean()),
                'max': float(masked.max()),
                'std': float(masked.std()),
                'sum': float(masked.sum()),
                'count': int(masked.count()),
                'median': float(np.ma.median(masked))}
                # 'mode': float(stats.mode(masked)[0][0])

            no_stats = {
                'min': -9999.,
                'mean': -9999.,
                'max': -9999.,
                'std': -9999.,
                'sum': -9999.,
                'count': -9999.,
                'median': -9999.}

            if np.nanmean(masked) == '--':
                statDict[feat.GetField(fldname)] = no_stats
            else:
                statDict[feat.GetField(fldname)] = feature_stats

    #clearing memory
        rvds = None
        mem_ds = None
    vds = None
    rds = None

    # print(statDict)
    
    ##OUTPUT
    #writes all stats to txt file (requires filenm)
    fields = ['staid', 'count', 'min', 'mean', 'median', 'max', 'std', 'sum']
    with open(filenm, 'wb') as f:
        w = csv.DictWriter(f, fields)
        w.writeheader()
        for key, val in statDict.items():
            # print(key,val)
            row = {'staid':key}
            row.update(val)
            w.writerow(row)
    f.close()


 
def zonal_category(ws, gdb, fc, fldname, ras, buffDist, cmap, filenm, nodata_value=None):
    #Open raster
    rds = gdal.Open(ras, GA_ReadOnly)
    assert(rds)
    rb = rds.GetRasterBand(1)
    rgt = rds.GetGeoTransform()
    rsize = (rds.RasterXSize, rds.RasterYSize)
    #Get raster cell size
   #this works in case cell size is not square
   #for square rasters len1 = len2
    dx1 = rgt[1]
    dy1 = rgt[2]
    dx2 = rgt[4]
    dy2 = rgt[5]
    len1 = np.sqrt(dx1**2 + dy1 **2)
    len2 = np.sqrt(dx2**2 + dy2 **2)
    ras_area = len1 * len2
    #if needed, angle of rotation (rotated rasters)
    theta = np.arctan(dy1/dx1)

    #Establish NoData value -- note: can be provided as input  
    if nodata_value:
        nodata_value = float(nodata_value)
        rb.SetNoDataValue(nodata_value)
    else:
        nodata_value = rb.GetNoDataValue()
    print('Raster NODATA Value: ', nodata_value)
    #Open feature class
    vds = ogr.Open(gdb, GA_ReadOnly)
    vlyr = vds.GetLayerByName(fc)
    
    #Create memory drivers to hold arrays
    mem_drv = ogr.GetDriverByName('Memory')
    driver = gdal.GetDriverByName('MEM')
 
    # Loop through vectors
    stats = []
    statDict = {}
    lenid = len(vlyr)

    for i, feat in enumerate(vlyr):
        fldid = feat.GetField(fldname)
        print('{} of {}, staid: {}'.format(i+1, lenid, fldid))
        #Buffer well points, using buffDist input (in meters)
        geom = feat.GetGeometryRef()
        buff = geom.Buffer(buffDist) 
        vec_area = buff.GetArea()
        # print('Ratio Buffer to Raster: ', vec_area/ras_area)

        src_offset = bbox_to_pixel_offsets(rgt, buff.GetEnvelope(),rsize)
        if src_offset[2] <= 0 or src_offset[3] <=0:
            #if point falls outside raster grid, include nodata as zone analysis
            statDict[feat.GetField(fldname)] = nodata_value
        else:    
            if vec_area/ras_area >= 4:
                zooms = 1
            else:
                zooms = int(np.sqrt(ras_area/(vec_area/4)))

            src_array = rb.ReadAsArray(*src_offset)
            # print(src_array)
            # plt.imshow(src_array)
            # plt.show()

            # calculate new geotransform of the feature subset
            new_gt = ((rgt[0] + (src_offset[0] * rgt[1])),
                rgt[1]/zooms,
                0.0,
                (rgt[3] + (src_offset[1] * rgt[5])),
                0.0,
                rgt[5]/zooms)
            
            # Create a temporary vector layer in memory
            mem_ds = mem_drv.CreateDataSource('out')
            mem_layer = mem_ds.CreateLayer('poly', None, ogr.wkbPolygon)
            mem_poly = ogr.Feature(mem_layer.GetLayerDefn())
            mem_poly.SetGeometryDirectly(buff)
            mem_layer.CreateFeature(mem_poly)
            
            # Rasterize the polygon
            rvds = driver.Create('', src_offset[2]*zooms, src_offset[3]*zooms, 1, gdal.GDT_Byte)
            rvds.SetGeoTransform(new_gt)
            gdal.RasterizeLayer(rvds, [1], mem_layer, burn_values=[1])
            rv_array = rvds.ReadAsArray()
            # print(rv_array)
            # plt.imshow(rv_array)
            # plt.show()

            #Resample the raster (only changes when zooms not 1)
            src_re = zoom(src_array, zooms, order = 0)
            # print(src_re)
            # plt.imshow(src_re)
            # plt.show()

            # Mask the source data array with our current feature
            # we take the logical_not to flip 0<->1 to get the correct mask effect
            # we also mask out nodata values explictly
            masked = np.ma.MaskedArray(src_re,
                mask=np.logical_or(src_re == nodata_value,
                    np.logical_not(rv_array)))
            # print(masked)
            # plt.imshow(masked)
            # plt.show()

            keys, counts = np.unique(masked.compressed(), return_counts=True)
            pixel_count = dict(zip([cmap[k] for k in keys],
                              [np.asscalar(c) for c in counts]))

            pixtot = float(sum(pixel_count.values()))

            for k, v in pixel_count.items():
                pixel_count[k] = v / pixtot * 100

            statDict[feat.GetField(fldname)] = pixel_count

    #clearing memory
        rvds = None
        mem_ds = None
 
    vds = None
    rds = None

    ##OUTPUT options
    #writes all stats to txt file (requires filenm)
    with open(filenm, 'w') as f:
        clss_all = cmap.values()
        # print(clss_all)
        fields = ','.join(['staid'] + clss_all)
        f.write('{}\n'.format(fields))
        for key, val in statDict.items():
            # print(key,val)
            line = []
            line.append(key)
            # print(key)
            # print(val)
            for v in clss_all:
                if v in val.keys():
                    line.append(str(val[v]))
                else:
                    line.append('0') #0 means there is none of that category in buffer
            line = ','.join(line)
            f.write('{}\n'.format(line))
    f.close()

