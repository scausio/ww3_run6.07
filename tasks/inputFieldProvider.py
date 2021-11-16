import xarray as xr
import os
import numpy as np
import glob
import pandas as pd
import subprocess
from datetime import datetime,timedelta,date
from glob import glob
from utils import dayBefore,dataFormatter,checkFunctionCompleted,seaoverland,getDaysBetweenDates, kelvin_to_celsius
from natsort import natsorted


def selectInputFile(conf,inputType,day):
    forecast = buildFilePath(conf, inputType, day, 'forecast')
    forecast = findLastForecast(forecast, day, inputType,conf)
    analysis = buildFilePath(conf, inputType, day, 'analysis')
    return forecastAnalysisPath(forecast, analysis)

def checkArraySize(air,sst):
    if air['T2M'].values.shape[0]>sst.values.shape[0]:
        t_out=copyDailyMean(sst.values,air['T2M'].values.shape)
    elif air['T2M'].values.shape[0]<sst.values.shape[0]: 
        t_out=sst.interp(time=[datetime.strptime(t,'%Y%m%d %H%M%S')for t in air['time']],method='nearest')
    else:
        t_out=sst
    return t_out

def copyDailyMean(dailymean,sizeTarget):
    timesInput=dailymean.shape[0]
    out=np.zeros(sizeTarget)*np.nan 
    print (out.shape)
    reps=int(sizeTarget[0]/timesInput)
    n=0
    for t in dailymean:
        for i in range(reps):
            out[n,:,:]=t
            n+=1
    return out

def cutField(conf, fieldName, lon, lat):
    tempFold = os.getcwd()
    cut = conf.cutArea

    if not os.path.exists(tempFold):
        os.makedirs(tempFold)
    latId = (lat >= cut.lat[0]) & (lat <= cut.lat[1])
    lonId = (lon >= cut.lon[0]) & (lon <= cut.lon[1])

    return lonId, latId

def findLastForecast(forecast,day,inputType,conf):
    try:
        glob(forecast)[0]
        return forecast
    except:
        n = 1
        while n <= 6:
            yesterday = (dataFormatter(day)- timedelta(days=n)).strftime('%Y%m%d')
            forecast = buildFilePath(conf, inputType, day, 'forecast',yesterday)
            print ('check forecast',forecast)
            try:
                glob(forecast)[0]
                return forecast
            except:
                n += 1

def buildFilePath(conf, dataset,day, analysisORforecast,yesterday=False):
    if dataset in ['waterVelocity','waterLevel','surfaceTemp']:
        fileconf = conf.copernicusFiles.parentHydro.fileConf
        filename = fillPathTemplates(os.path.join(fileconf[analysisORforecast].base,fileconf[analysisORforecast].nameTemplate), fileconf, day, yesterday)
    else:
        fileconf = conf.copernicusFiles[dataset].fileConf
        filename = fillPathTemplates(os.path.join(fileconf[analysisORforecast].base,fileconf[analysisORforecast].nameTemplate), fileconf, day, yesterday)
    return filename

def forecastAnalysisPath(forecast, analysis):
    print ('forecast to check:',forecast)
    print ('analysis to check:',analysis)
    try:# analysis
        path=natsorted(glob(analysis))
        path[0]
    except:  # forecast
        path=natsorted(glob(forecast))
        path[0]

    print(path, ':SELECTED')
    return path

def saveGridSpecs(lon1, lon2, lat1, lat2, lonSize, latSize, v):
    lines = [lon1, lon2, lonSize, lat1, lat2, latSize]
    with open('{}_gridSpecs.txt'.format(v), 'w') as out:
        out.write(" '{}' 'LL' T T\n".format(v))
        for l in lines:
            out.write(' {} '.format(l))
        out.close()


def fillPathTemplates(template, fileConf, startdate,yesterday):
    d = {'date': str(startdate),
         'year': str(startdate)[:4],
         'month': str(startdate)[4:6],
         'day': str(startdate)[6:8],
         'yesterday': dayBefore(str(startdate)) if not yesterday else yesterday,
         'freq': fileConf.freq,
         'producer': fileConf.producer,
         'parameter': fileConf.parameter,
         'config': fileConf.config,
         'region': fileConf.region,
         'version': fileConf.version
         }
    print (template)
    return template.format(d=d)


def getHydroData(path, conf, fieldType):
    #field type can be: waterVelocity,waterLevel,surfaceTemp
    if os.path.exists('HYDRO.nc'):
        print ('hydro already in the directory')

    fo = xr.open_mfdataset(path,combine='by_coords')
    lat,lon = fo[conf.copernicusFiles.parentHydro.lat][:], fo[conf.copernicusFiles.parentHydro.lon][:]

    # get data inside box
    lonId, latId = cutField(conf, 'oceanMsk', lon, lat)
    fo_sub = fo.isel(**{conf.copernicusFiles.parentHydro.lat:np.where(latId)[0], conf.copernicusFiles.parentHydro.lon:np.where(lonId)[0], conf.copernicusFiles.parentHydro.depth:0})
    #fo[conf.copernicusFiles.parentHydro.time] = fo[conf.copernicusFiles.parentHydro.time].dt.floor('H')  # round hours

    if fieldType=='waterVelocity':
        print ('getting currents')
        u = fo_sub[conf.copernicusFiles.parentHydro.waterVelocity.u]
        u.values[u.values==0]=np.nan
        u.values = np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)), 20) for f in u.values])
        print ('u current sol', u.values.shape)
        v = fo_sub[conf.copernicusFiles.parentHydro.waterVelocity.v]
        v.values[v.values==0]=np.nan
        v.values = np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)), 20) for f in v.values])
        print ('v current sol', v.values.shape)
        fo_sub[conf.copernicusFiles.parentHydro.waterVelocity.u].values = u
        fo_sub[conf.copernicusFiles.parentHydro.waterVelocity.v].values = v
        try:
            fo_sub = fo_sub.compute()
        except:
            pass
        fo_sub.to_netcdf('CUR.nc')

    if fieldType == 'waterLevel':
        lev= fo_sub[conf.copernicusFiles.parentHydro.waterLevel.ssh]
        print ('lev', lev.values.shape)
        lev.values[lev.values==0]=np.nan
        lev.values = np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)), 20) for f in lev.values])
        fo_sub = fo.isel(lat=latId, lon=lonId)
        fo_sub[conf.copernicusFiles.parentHydro.waterLevel.ssh].values = lev
        try:
            fo_sub = fo_sub.compute()
        except:
            pass
        fo_sub.to_netcdf('LEV.nc')

    if fieldType == 'surfaceTemp':
        sst = fo_sub[conf.copernicusFiles.parentHydro.surfaceTemp.sst]
        sst.values[sst.values==0]=np.nan
        sst.values = np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)), 20) for f in sst.values])
        print ('sst current sol', sst.values.shape)
        fo_sub[conf.copernicusFiles.parentHydro.surfaceTemp.sst].values = sst
        try:
            fo_sub = fo_sub.compute()
        except:
            pass

    return fo_sub


def getMeteoData(path, conf):
    #fo = xr.open_mfdataset(path,combine='by_coords')
    fo = xr.open_mfdataset(path,combine='by_coords')#.sel({conf.copernicusFiles.meteoData.lon:slice(conf.cutArea.lon[0],conf.cutArea.lon[0]), conf.copernicusFiles.meteoData.lat:slice(conf.cutArea.lat[0],conf.cutArea.lat[1])})
    lon, lat = fo[conf.copernicusFiles.meteoData.lon][:], fo[conf.copernicusFiles.meteoData.lat][:]
   
    lonId, latId = cutField(conf, 'meteoMsk', lon, lat)
    #fo_sub = fo.isel(**{conf.copernicusFiles.parentHydro.lat:np.where(latId)[0], conf.copernicusFiles.parentHydro.lon:np.where(lonId)[0]})
    land=fo['LSM'][:, latId, lonId].compute()
    #land=fo_sub['LSM'][:]
    msk_land=land.data>0.5
    msk_land=np.logical_not(msk_land).astype(int)
    
    #fo[conf.copernicusFiles.meteoData.time] = fo[conf.copernicusFiles.meteoData.time].dt.floor('H')  # round hours
    u=fo[conf.copernicusFiles.meteoData.u][:, latId, lonId].compute()
    #u=fo_sub[conf.copernicusFiles.meteoData.u][:].compute()
    u_data=u.values
    u_data=u_data*msk_land
    u_data[u_data==0]=np.nan
    print ('u', u.shape)
    u_data=np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)),20)for f in u_data])
    print ('sol u completed')

    v=fo[conf.copernicusFiles.meteoData.v][:, latId, lonId].compute()
    #v=fo_sub[conf.copernicusFiles.meteoData.v][:].compute()
    v_data=v.values
    v_data=v_data*msk_land
    v_data[v_data==0]=np.nan
    v_data=np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)),20)for f in v_data])
    print ('sol v completed')

    t2m=fo[conf.copernicusFiles.meteoData.T2M][:, latId, lonId].compute()
    #t2m=fo_sub[conf.copernicusFiles.meteoData.T2M][:].compute()
    t2m_data=t2m.values
    t2m_data=t2m_data*msk_land
    t2m_data[t2m_data==0]=np.nan
    t2m_data=np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)),20)for f in t2m_data])
    print ('sol t2m completed')

    fo_sub=fo.isel(lat= latId,lon= lonId)
    fo_sub[conf.copernicusFiles.meteoData.v].values  = v_data #fix
    fo_sub[conf.copernicusFiles.meteoData.u].values  = u_data #fix
    fo_sub[conf.copernicusFiles.meteoData.T2M].values = t2m_data

    fo_sub=fo_sub.compute()
    fo_sub=fo_sub.sortby(conf.copernicusFiles.meteoData.lat, ascending=True)
    fo_sub.to_netcdf('WND.nc')

    return fo_sub

class WNSgetter():
    def __init__(self, conf, date,duration, fillValue):
        self.workPath = os.getcwd()
        self.conf = conf
        self.date = date
        self.fillValue = fillValue
        self.duration=int(duration)
        print('getting WNS')

    def manageSST(self):

        date = dataFormatter(self.date)
        day_before=(date - timedelta(days=1)).strftime('%Y%m%d')
        days=getDaysBetweenDates(day_before,self.duration+2)

        sstFiles=[]
        for day in days:
            ssts=selectInputFile(self.conf, 'surfaceTemp', day)
            if isinstance(ssts, list):
                for sst in ssts:
                    if sst not in sstFiles:
                        sstFiles.append(sst)
            else:
                sstFiles.append(ssts)



        #sstFiles = [selectInputFile(self.conf, 'surfaceTemp', day) for day in days]
        SST= getHydroData(sstFiles, self.conf, 'surfaceTemp')
        return SST

    def manageWind(self):
        date = dataFormatter(self.date)
        day_before=(date - timedelta(days=1)).strftime('%Y%m%d')

        days=getDaysBetweenDates(day_before,self.duration+2)


        windFiles = []
        for day in days:
            wnds = selectInputFile(self.conf, 'meteoData', day)
            if isinstance(wnds, list):
                for wnd in wnds:
                    if wnd not in windFiles:
                        windFiles.append(wnd)
            else:
                windFiles.append(wnds)

        #wndInputFiles=[selectInputFile(self.conf,'meteoData',day) for day in days]
        WND= getMeteoData(windFiles, self.conf)
        return WND

    def getDeltaT(self, wind, hydro,conf):
        print ('computing DT')

        #wind['T2M'].values[wind['T2M'].values==self.fillValue]=np.nan
        
        sst=hydro[conf.copernicusFiles.parentHydro.surfaceTemp.sst] 
        sst.values[sst.values==0]=np.nan
        sst.values= np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)), 20) for f in sst])

        regridded_T=sst.interp(**{conf.copernicusFiles.parentHydro.lat:wind[conf.copernicusFiles.meteoData.lat].values,
                                                                          conf.copernicusFiles.parentHydro.lon:wind[conf.copernicusFiles.meteoData.lon].values})
        sea_T=regridded_T.interp(method='nearest',**{conf.copernicusFiles.parentHydro.time: wind[conf.copernicusFiles.meteoData.time]})
        #conf.copernicusFiles.parentHydro.time: conf.copernicusFiles.meteoData.time}
        #nearest_T.values[nearest_T.values==0]=np.nan
        sea_T= np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)), 20) for f in sea_T])
        air_T= wind[conf.copernicusFiles.meteoData.T2M] 
        air_T.values[air_T.values==self.fillValue]=np.nan
        if air_T.attrs['units'] in ['K', 'Kelvin', 'k', 'degK']:
            air_T=kelvin_to_celsius(air_T)
        dT = air_T.values- sea_T
        wind["DT"]=((conf.copernicusFiles.meteoData.time, conf.copernicusFiles.meteoData.lat, conf.copernicusFiles.meteoData.lon),dT)
        wind["DT"].attrs["units"]="Celsius"
        wind=wind.drop(conf.copernicusFiles.meteoData.T2M)
        wind.to_netcdf('WNS.nc')

    def getWNS(self):
        WND = self.manageWind()  # ,T2M,latW,lonW
        SST = self.manageSST()  # ,latT,lonT # u can force computation with force=True
        self.getDeltaT(WND, SST,self.conf)


class WNDgetter():
    def __init__(self, conf, date,duration, fillValue):
        self.workPath = os.getcwd()
        self.conf = conf
        self.date = date
        self.fillValue = fillValue
        self.duration=int(duration)
        print('getting WND')

    def manageWind(self):

        days=getDaysBetweenDates(self.date,self.duration+1)

        windFiles = []
        for day in days:
            wnds = selectInputFile(self.conf, 'meteoData', day)
            if isinstance(wnds, list):
                for wnd in wnds:
                    if wnd not in windFiles:
                        windFiles.append(wnd)
            else:
                windFiles.append(wnds)

        #wndInputFiles=[selectInputFile(self.conf,'meteoData',day)for day in days]
        WND=getMeteoData(windFiles, self.conf)
        return WND

    def getWND(self):
        WND = self.manageWind()
        print('... done ...')


class CURgetter():

    def __init__(self, conf, date,duration, fillValue):
        self.workPath = os.getcwd()
        self.conf = conf
        self.date = date
        self.fillValue = fillValue
        self.duration=int(duration)
        print('getting Currents')


    def manageCUR(self):
        days=getDaysBetweenDates(self.date,self.duration+1)
        curFiles=[]
        for day in days:
            curs = selectInputFile(self.conf, 'waterVelocity', day)
            if isinstance(curs, list):
                for cur in curs:
                    if cur not in curFiles:
                        curFiles.append(cur)
            else:
                curFiles.append(curs)

        #curFiles = [selectInputFile(self.conf, 'waterVelocity', day) for day in days]
        CUR= getHydroData(curFiles, self.conf, 'waterVelocity')
        return CUR

    def getCUR(self):
        CUR = self.manageCUR()  # ,latT,lonT # u can force computation with force=True
        print('... done ...')


class LEVgetter():
    def __init__(self,  conf, date,duration, fillValue):
        self.workPath = os.getcwd()
        self.conf = conf
        self.date = date
        self.fillValue = fillValue
        self.duration=int(duration)
        print('getting Sea Level')


    def manageLEV(self):
        days=getDaysBetweenDates(self.date,self.duration+1)

        levFiles = []
        for day in days:
            levs = selectInputFile(self.conf, 'waterLevel', day)
            if isinstance(levs, list):
                for lev in levs:
                    if lev not in levFiles:
                        levFiles.append(lev)
            else:
                levFiles.append(levs)

        #levFiles = [selectInputFile(self.conf, 'waterLevel', day) for day in days]
        LEV, time, lat, lon = getHydroData(levFiles, self.conf, 'waterLevel')
        return {'LEV': LEV, 'lat': lat, 'lon': lon, 'time': time}

    def getLEV(self):
        outname = 'LEV.data'
        LEV = self.manageLEV()
        print('... done ...')
