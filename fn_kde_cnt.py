#JTVD - 2020

#import libraries
import pandas as pd
import geopandas as gpd
import numpy as np
import shapely.wkt
import sys

#classes, method
from osgeo import ogr
from hull.concavehull import ConcaveHull
from sklearn.cluster import dbscan
from shapely.geometry import Polygon

#capture stderr
sys.stderr = object

#concave points
def to_concave_points(df,coord):
    groups = np.unique(df['group'].tolist())
    contours = []
    for g in groups:
        gdf = df.loc[df['group'] == g]
        coords = [[int(x[0]),int(x[1])] for x in (list(zip(gdf.x,gdf.y)))]
        concave = ConcaveHull(coords)
        concave_array = concave.calculate()
        if not concave_array is None:
            contours.append(concave_array)
    #return contours list
    return contours

#concave polygons
def to_concave_polygons(df,lvl,grid):

    #select levels
    val = df[(df['kde'] == lvl)]

    #not empty
    if not val.empty:

        #reconstruct grid
        grd = pd.merge(grid,val,on='gid',how='inner')

        #cluster
        coord = [[int(x[1]),int(x[0])] for x in (list(zip(grd.x,grd.y)))]
        cs,lbls = dbscan(coord,eps=2000)
        grd['group'] = lbls

        #group to concave points
        contourpnt = to_concave_points(df=grd[(grd['group'] >= 0)],coord=coord)

        #select
        contourset = gpd.GeoSeries([Polygon(contour) for contour in contourpnt if len(contour) >= 10])
        contourset = gpd.GeoDataFrame({'geometry': contourset})
        contourset.crs = 'epsg:27700'

        #pretty
        contourset['geometry'] = contourset.geometry.buffer(10000,join_style=1).buffer(-10000,join_style=1)

        #meta
        contourset['level'] = lvl
        contourset['idx'] = contourset.index
        return(contourset)

    #empty
    else:
        return(val)

#read from stdin
input = pd.read_csv(sys.stdin,sep=';',names=['surname','year','freq','bw','idx','kde'])

#read from disk
gridc = pd.read_csv('data/xy.csv').sort_values(by='gid')
gbshp = gpd.read_file('data/gb.shp')
gbshp.crs = 'epsg:27700'
gbwkt = ogr.CreateGeometryFromWkt(str(gbshp['geometry'][0]))

#kde df
idx = [int(x) for x in str(input['idx'].to_numpy())[2:-2].split(',')]
kde = [int(x) for x in str(input['kde'].to_numpy())[2:-2].split(',')]
kdf = pd.DataFrame({'gid':idx,'kde':kde})

#levels
lvl1 = to_concave_polygons(df=kdf,lvl=1,grid=gridc)
lvl2 = to_concave_polygons(df=kdf,lvl=2,grid=gridc)
lvl3 = to_concave_polygons(df=kdf,lvl=3,grid=gridc)

#overlay1
try:
    mrg1 = gpd.overlay(lvl1,lvl2,how='difference')
except AttributeError:
    if not lvl1.empty and lvl2.empty:
        mrg1 = lvl1
    elif not lvl2.empty and lvl1.empty:
        mrg1 = lvl2
    else:
        mrg1 = pd.DataFrame()

#overlay2
try:
    mrg2 = gpd.overlay(lvl2,lvl3,how='difference')
except AttributeError:
    if not lvl2.empty and lvl3.empty:
        mrg2 = lvl2
    elif not lvl3.empty and lvl2.empty:
        mrg2 = lvl3
    else:
        mrg2 = pd.DataFrame()

#combine
comb1 = mrg1.append(mrg2)
comb2 = comb1.append(lvl3)

#output
try:

    #harmonise
    comb2.drop_duplicates(subset=['idx','level'],inplace=True,keep='first')
    comb2.reset_index(inplace=True)
    kdeshp = comb2.reset_index()

    #gb intersect through osgeos
    for index,row in kdeshp.iterrows():
        cnwkt = ogr.CreateGeometryFromWkt(str(row['geometry']))
        kdeshp.at[index,'intersect'] = cnwkt.Intersection(gbwkt).ExportToWkt()

    #cleanup
    kdeshp['geometry'] = kdeshp['intersect'].map(shapely.wkt.loads)
    kdeshp.drop(kdeshp.columns[[0,1,4,5]],axis=1,inplace=True)
    kdeshp.crs = 'epsg:27700'

    #reproject
    kdeprj = kdeshp.to_crs('epsg:4326')

    #output
    kdejson = str(kdeprj.to_json())
    print(input.iloc[0]['surname']+';'+str(input.iloc[0]['year'])+';'+str(input.iloc[0]['freq'])+';'+str(input.iloc[0]['bw'])+';'+kdejson)

except:

   #output
   print(input.iloc[0]['surname']+';'+str(input.iloc[0]['year'])+';'+str(input.iloc[0]['freq'])+';'+str(input.iloc[0]['bw'])+';NULL')
