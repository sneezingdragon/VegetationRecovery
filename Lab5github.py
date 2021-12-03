import rasterio
import numpy as np
import pandas as pd
from math import pi
from scipy import ndimage
import glob
import scipy
import matplotlib.pyplot as plt
fire = rasterio.open(r'C:\Users\sneez\OneDrive\Documents\lab5\data\fire_perimeter.tif')
dem = rasterio.open(r'C:\Users\sneez\OneDrive\Documents\lab5\data\bigElk_dem.tif')
demdata = dem.read(1)
cellsize = 30
def slopeAspect(dem, cs):
    kernel = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]])
    dzdx = ndimage.convolve(dem, kernel, mode='mirror') / (8 * cs)
    dzdy = ndimage.convolve(dem, kernel.T, mode='mirror') / (8 * cs)
    slp = np.arctan((dzdx ** 2 + dzdy ** 2) ** 0.5) * 180 / pi
    ang = np.arctan2(-dzdy, dzdx) * 180 / pi
    aspect = np.where(ang > 90, 450 - ang, 90 - ang)
    return slp, aspect

def reclassAspect(npArray):

    return np.where((npArray > 22.5) & (npArray <= 67.5), 2,
    np.where((npArray > 67.5) & (npArray <= 112.5), 3,
    np.where((npArray > 112.5) & (npArray <= 157.5), 4,
    np.where((npArray > 157.5) & (npArray <= 202.5), 5,
    np.where((npArray > 202.5) & (npArray <= 247.5), 6,
    np.where((npArray > 247.5) & (npArray <= 292.5), 7,
    np.where((npArray > 292.5) & (npArray <= 337.5), 8, 1)))))))

def reclassByHisto(npArray, bins):
    histo = np.histogram(npArray, bins)[1]
    rClss = np.zeros_like(npArray)
    for i in range(bins):
        rClss = np.where((npArray >= histo[i]) & (npArray <= histo[i + 1]),
                         i + 1, rClss)
    return rClss 
slopeasp = slopeAspect(demdata, cellsize)
slope, asp= slopeasp
cardinalaspect = reclassAspect(asp)
slopereclass = reclassByHisto(slope, 10)
def NDVI(year):
    holdmyspot= []
    for files in glob.glob(r'C:\Users\sneez\OneDrive\Documents\lab5\data\L5_big_elk\*.tif'):
        if files[-11:-7] == year:
            rasters = rasterio.open(files)
            
            if files[-6:-4] == 'B3':
                B3= rasters.read(1)
                
            if files[-6:-4] == 'B4':
                B4 = rasters.read(1)
                
                NDVIyear = np.divide((np.subtract((B4), (B3))), (np.add((B4), (B3))))
                return(NDVIyear)
year = (['2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011'])
for x in year:
    NDVI(x)
    print('we did it!')
fireraster = fire.read(1)
FireyBool = np.where(fireraster == 2, 1, 0)
year = np.array(['2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011'], dtype= 'float64')
def RecRat(year):
    temp_arr= np.zeros_like(NDVI(year))
    healthyforest = np.multiply(NDVI(year), FireyBool)
    nonzero = (healthyforest[healthyforest != 0])
    meanofyear = np.mean(nonzero)
    RRyears= np.divide(NDVI(year), meanofyear)
    flatstanley = RRyears.flatten()
    grocerylist= list(flatstanley)
    return grocerylist
RRvalues = np.column_stack([RecRat('2002'), RecRat('2003'), RecRat('2004'), RecRat('2005'), RecRat('2006') 
                            , RecRat('2007'), RecRat('2008'), RecRat('2009'), RecRat('2010'), RecRat('2011')])
AllRRvalues= np.transpose(RRvalues)
def OnlyBurned(year):
    testyear= (RecRat(year)[FireyBool ==0])
    return testyear
RRline =np.polyfit(year, AllRRvalues, 1)
onlyslope = np.delete(RRline, 1, 0)
backtoblack = np.reshape(onlyslope, (280, 459))
healthy = backtoblack[FireyBool ==0]
print(('The mean coefficient of recovery is'), healthy.mean())
def RRMean(year):
    pirates = []
    meanrr = sum(RecRat(year))/len(RecRat(year))
    pirates.append(meanrr)
    return pirates
year = (['2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011'])
for x in year:
    print(RRMean(x))
    def ZonalStats(RecovRatio, SlopeorAss, csvname):
    specialzones = np.unique(SlopeorAss)
    dic={'zones':[], 'mean':[], 'max':[], 'min':[], 'sd':[], 'count':[]}
    x=1
    Recover= np.asarray(RecovRatio)
    Recovery=np.reshape(Recover, (280, 459))
    for zones in list(specialzones):
        dic['zones'].append(x)
        oneatatime = np.where(SlopeorAss==x,1,np.nan)
        dic['mean'].append(np.nanmean(oneatatime * RecovRatio))
        dic['max'].append(np.nanmean(oneatatime * RecovRatio))
        dic['min'].append(np.nanmean(oneatatime * RecovRatio))
        dic['sd'].append(np.nanstd(oneatatime * RecovRatio))
        dic['count'].append(np.nansum(oneatatime))
        x=x+1
    df = pd.DataFrame(dic)
    csvname=df.to_csv(csvname)
ZonalStats((backtoblack), slopereclass, 'SlopeStatsAll.csv')
ZonalStats((backtoblack), cardinalaspect, 'AspectStatsAll.csv')
reverse = np.where(fireraster == 2, 0, 1)
yaburnt= np.where(reverse ==1, backtoblack, 0)
yaburnt[140, 220]
with rasterio.open(r'C:\Users\sneez\OneDrive\Documents\lab5\data\Finito.tif', 'w',
                        driver = "GTiff",
                        height = yaburnt.shape[0],
                        width = yaburnt.shape[1],
                        count = 1,
                        dtype = 'float32',
                        crs = fire.crs,
                        transform = fire.transform,
                        nodata = 0,
                    )as out_raster:
                    out_raster.write(yaburnt, 1)
        
print("Areas with the worst coefficient of recovery were on North facing slopes, while vegetation on South facing \n"
"slopes had a better recovery rate over the years analyzed. Overall recovery was better on areas with less slope, at \n" 
"low elevations, probably due to erosion and runoff. Precipitation levels can effect regrowth, as well as absorption of \n"
"water in soil beds. \n")