# -*- coding: utf-8 -*-
'''
Name:        a.fpcreator.FloodPlainDelineation
'''

import os
import sys
import json

#raster statistic libraries
import rasterio
from   rasterio.mask    import mask
from   shapely.geometry import mapping

#hidden import for pyinstaller
#import pkg_resources.py2_warn
from pyproj import _datadir, datadir
from fiona import _shim, schema

#standard module
import pandas           as     pd
import numpy            as     np
import geopandas        as     gpd

#shapely aus (Shapely, Geos, und Fiona erstellt)
from   shapely.geometry import Polygon,shape,LinearRing,Point,MultiPolygon
from   shapely.ops      import unary_union
from   rasterio.crs     import CRS

# gok = r"y:\ZZ_Entwicklung\FP_creator\mxd\Input_Hydraulik\shp\KANTEN3D_04TestHW2014_g1.shp"
# wsp = r"y:\ZZ_Entwicklung\FP_creator\mxd\Input_Hydraulik\shp\KANTEN3D_04TestHW2014_w1.shp"

# df_gok = gpd.read_file(gok)
# df_wsp = gpd.read_file(wsp)

wsp_max = r"y:\Hyd2D_Ooser_Landgr_DB\M01_Hydro_AS-2D\02_Ist-Zustand\HQ1000_48h\wspl_max.dat"
netz    = r"y:\Hyd2D_Ooser_Landgr_DB\M01_Hydro_AS-2D\02_Ist-Zustand\HQ1000_48h\hydro_as-2d.2dm"

wmax = np.loadtxt(wsp_max,skiprows=4,dtype=np.float32)
node = np.arange(1,len(wmax)+1,1)

idx   = np.where(wmax!=0)
aktiv = node[idx]

elem = np.loadtxt(netz,usecols=0,dtype=np.str,  skiprows=3,max_rows=814520)
eid  = np.loadtxt(netz,usecols=1,dtype=np.int32,skiprows=3,max_rows=814520)
e1   = np.loadtxt(netz,usecols=2,dtype=np.int32,skiprows=3,max_rows=814520)
e2   = np.loadtxt(netz,usecols=3,dtype=np.int32,skiprows=3,max_rows=814520)
e3   = np.loadtxt(netz,usecols=4,dtype=np.int32,skiprows=3,max_rows=814520)
e4   = np.loadtxt(netz,usecols=5,dtype=np.int32,skiprows=3,max_rows=814520)

#triangles
tri  = np.where(elem == 'E3T')[0]
triangles = np.dstack((e1,e2,e3))[0][tri]

#rectangles
rec  = np.where(elem == 'E4Q')[0]
rectangles = np.dstack((e1,e2,e3,e4))[0][rec]

#read nodes
nodes  = np.loadtxt(netz,usecols=1,dtype=np.int32, skiprows=814523,max_rows=439931)
coords = np.loadtxt(netz,usecols=[2,3,4],dtype=np.float32, skiprows=814523,max_rows=439931)

pols = [*[Polygon([coords[i-1] for i in [*t,t[0]]]) for t in triangles],
        *[Polygon([coords[i-1] for i in [*r,r[0]]]) for r in rectangles]]

df  = pd.DataFrame({'geometry':pols})
gdf = gpd.GeoDataFrame(df)


xmin,xmax = coords[:,0].min(),coords[:,0].max()
ymin,ymax = coords[:,1].min(),coords[:,1].max()

x_grid = np.arange(int(xmin),int(xmax+1.5),0.5)
y_grid = np.arange(int(ymin),int(ymax+1.5),0.5)

x,y = np.meshgrid(x_grid,y_grid)

# def nodes2coords(t1,t2,t3):
#     return Polygon([coords[i-1] for i in [t1,t2,t3,t1]])

# p = np.vectorize(nodes2coords)(triangles[:,0],triangles[:,1],triangles[:,2]).tolist()
# netz3T = MultiPolygon(p)
