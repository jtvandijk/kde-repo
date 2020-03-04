#JTVD - 2020

#import libraries
import pandas as pd
import geopandas as gpd
import numpy as np
import sys

#classes, method
from sklearn.cluster import dbscan
from shapely.geometry import Point,Polygon,MultiPolygon
from fiona.crs import from_epsg

#capture stderr
sys.stderr = object

#coordinate lookup
def create_lookup(coord):
    record = dict()
    sorted_data = sorted(coord)
    for x,y in sorted_data:
        if x in record.keys():
            record[x][y] = None
        else:
            record[x] = dict()
            record[x][y] = None
    return record

#concave points
def to_concave_points(df,coord):
    groups = np.unique(df['group'].tolist())
    contours = []
    for g in groups:
        boundary = []
        gdf = df.loc[df['group'] == g]
        X = np.unique(gdf['x'].tolist())
        coordx = [[int(x[0]),int(x[1])] for x in (list(zip(gdf.x,gdf.y)))]
        Xlook = create_lookup(coord=coordx)
        for x in X:
            if x in Xlook.keys():
                min_y = min(Xlook[x].keys())
                max_y = max(Xlook[x].keys())
                if not min_y == max_y:
                    boundary.append((int(x),min_y))
        for x in np.flipud(X):
            if x in Xlook.keys():
                max_y = max(Xlook[x].keys())
                min_y = min(Xlook[x].keys())
                if not max_y == min_y:
                    boundary.append((int(x),max_y))
        contours.append(boundary)
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
        contourset.crs = from_epsg(27700)

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
gbshp.crs = from_epsg(27700)

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

    #gb clip
    kdeshp = gpd.overlay(comb2,gbshp,how='intersection')

    #re-project
    kdeshp['geometry'] = kdeshp['geometry'].to_crs(epsg=4326)
    kdeshp.drop(kdeshp.columns[[1,2]],axis=1,inplace=True)

    #output
    kdejson = str(kdeshp.to_json())
    print(input.iloc[0]['surname']+';'+str(input.iloc[0]['year'])+';'+str(input.iloc[0]['freq'])+';'+str(input.iloc[0]['bw'])+';'+kdejson)

except:
    print(input.iloc[0]['surname']+';'+str(input.iloc[0]['year'])+';'+str(input.iloc[0]['freq'])+';'+str(input.iloc[0]['bw'])+';NULL')
