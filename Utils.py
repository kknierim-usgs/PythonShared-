from __future__ import print_function
__author__ = 'kknierim'

#import packages
import os,sys
import pandas as pd

#import arcpy
import arcpy
from arcpy import env

#General functions for converting between text files and ArcMap feature classes 

#convert feature class to textfile
def FCtoTXT(fc, txtloc, outtxt):
    env.workspace = txtloc
    arcpy.env.overwriteOutput = True
    fc_view = arcpy.MakeTableView_management(fc,'fc_view')
    arcpy.TableToTable_conversion(fc_view,txtloc,outtxt)
    del fc_view
    print("text file {} created".format(outtxt))

#convert text file to feature class
def TXTtoXYFC(ws, txt, sep, latfld, longfld, spatval, txtflds, fcloc, fcname):
    env.workspace = ws
    arcpy.env.overwriteOutput = True
    strDict = {}
    for f in txtflds:
        strDict[f] = str
    if sep == ',':
        data = pd.read_table(txt, sep = ',',dtype = strDict)
    else:
        data = pd.read_table(txt, dtype = strDict)
    print(data.head(5))
    y_coord = latfld
    x_coord = longfld
    spat = arcpy.SpatialReference(spatval)
    arcpy.MakeXYEventLayer_management(txt, x_coord, y_coord, 'outlyr', spat)
    outlyr = 'outlyr'
    print(arcpy.GetCount_management(outlyr))
    fms = arcpy.FieldMappings()
    fms.addTable(outlyr)
    for f in sitenoflds:
        fidx = data.columns.get_loc(f)
        fmap = fms.getFieldMap(fidx)
        ofield = fmap.outputField
        ofield.type = 'String'
        ofield.length = 30
        fmap.outputField = ofield
        fms.replaceFieldMap(fidx,fmap)
    print('field mapping complete')
    arcpy.FeatureClassToFeatureClass_conversion(outlyr, fcloc, fcname, "#", fms)
    print('feature class created, {}'.format(fcname))

#convert polygon to raster
def PolytoRas(ws, fc, fld, cellsize, outras):
    env.workspace = ws
    arcpy.env.overwriteOutput = True
    arcpy.PolygonToRaster_conversion(fc, fld, outras, 'MAXIMUM_COMBINED_AREA', '#', cellsize)
    print('feature class converted to raster')

#dictionary of useful spatial projections and WKIDs
def SpatVals():
    spatDict = {}
    spatDict['USGSAlbers'] = 102039 #USGS Albers
    spatDict['NAD83'] = 6318 #Geographic nad_83 (2011)
    return(spatDict)






