import os
import xarray as xr
import numpy as np
from utils import getConfigurationByID,getFreqs, getDirs, degToRad, duplicatesId, seaoverland
from toSpectra import JONSWAP
from utils import stringAllocation,myMFdataset, reverse90,fixBCdir
import netCDF4 as nc
from datetime import datetime


'''


* get header 1 time
    File ID in quotes, number of frequencies, directions and points.
    grid name in quotes (for unformatted file C*21,3I,C*30).
    Bin frequencies in Hz for all bins.
    Bin directions in radians for all bins (Oceanographic conv.).


- open parent model
- getNearest nodes
- for each node
    for each time step
        get jonswap
        get d, U10 and direction
        
        save txt with:
            - Time in yyyymmdd hhmmss format
            Point name , lat, lon, d, U10 and direction , current speed and direction 

'''


def getBCdataset(times,stat,freq,dir,outname,computed=False):

    ds_bc = nc.Dataset(outname, 'w', format='NETCDF4')
    ds_time = ds_bc.createDimension('time', len(times))
    ds_stat = ds_bc.createDimension('station', len([stat]))
    ds_freq=ds_bc.createDimension('frequency', len(freq))
    ds_dir=ds_bc.createDimension('direction', len(dir))
    string16 = ds_bc.createDimension('string16', 16)

    timess = ds_bc.createVariable('time', 'f8', ('time',))
    stats = ds_bc.createVariable('station', 'i4', ('station',))
    dirs = ds_bc.createVariable('direction', 'f4', ('direction',))
    freqs = ds_bc.createVariable('frequency', 'f4', ('frequency',))
    lats = ds_bc.createVariable('latitude', 'f4', ('time','station',))
    lons = ds_bc.createVariable('longitude', 'f4', ('time','station',))
    dpts = ds_bc.createVariable('depth', 'f4', ('time','station',))
    wss = ds_bc.createVariable('u10m', 'f4', ('time','station',))
    wdirs = ds_bc.createVariable('udir', 'f4', ('time','station',))
    curs = ds_bc.createVariable('curr', 'f4', ('time','station',))
    curdirs = ds_bc.createVariable('currdir', 'f4', ('time','station',))
    efth = ds_bc.createVariable('efth', 'f8', ('time', 'station','frequency', 'direction',))

    sn = ds_bc.createVariable('station_name', 'str', ('station','string16',))
    sn.long_name = "station name"
    sn.content = "XW"
    sn.associates = "station string16"

    st = ds_bc.createVariable('string16', 'i4', ('string16',))
    st.long_name = "station_name number of characters"
    st.axis = "W";

    dates = [(datetime.strptime(str(t)[:19], '%Y-%m-%dT%H:%M:%S') - datetime.fromtimestamp(0)).total_seconds()   for t in times.values]

    timess[:] = dates
    timess.units = "seconds since 1970-01-01 00:00:00.0"
    timess.calendar = "gregorian"
    timess.long_name = "time";
    timess.standard_name = "time";
    timess.axis = "T";

    if computed==True:
        dirs[:]=dir
    else:
        #in case of wam spectra bc
        dirs[:]=reverse90(dir)

    freqs[:]=freq
    stats[:]=stat
    st[:]=str(stat)
    curs[:]=np.zeros(len(times))
    curdirs[:]=np.zeros(len(times))

    stats.long_name= "station id"
    stats.axis= "X"

    lats.units="degree_north"
    lats.long_name="latitude"
    lats.standard_name = "latitude"
    lats.valid_min = -90.
    lats.valid_max = 90.
    lats.content = "TY"
    lats.associates = "time"

    lons.units="degree_east"
    lons.long_name="longitude"
    lons.standard_name = "longitude"
    lons.valid_min = -180.
    lons.valid_max = 180.
    lons.content = "TX"
    lons.associates = "time"

    freqs.units ="s-1"
    freqs.long_name="frequency of center band"
    freqs.standard_name ="sea_surface_wave_frequency"
    freqs.globwave_name="frequency"
    freqs.valid_min=0.
    freqs.valid_max=10.
    freqs.axis="Y"

    dirs.units ="degree"
    dirs.long_name="sea surface wave to direction"
    dirs.standard_name ="sea surface wave to direction"
    dirs.globwave_name="direction"
    dirs.valid_min=0.
    dirs.valid_max=360.
    dirs.axis="Z"
    #ds_bc['direction'].attrs['_FillValue'] = 9.96921e+36

    efth.units ="m2 s rad-1"
    efth.long_name="sea surface wave directional variance spectral density"
    efth.standard_name ="sea surface wave directional variance spectral density"
    efth.globwave_name="directional_variance_spectral_density"
    efth.valid_min=0.
    efth.valid_max=1.e+20
    efth.scale_factor=1.
    efth.add_offset = 0.
    efth.content = ""
    efth.associates = "time frequency direction"
    #ds_bc['efth'].attrs['_FillValue']=9.96921e+36
    efth.axis="TXYZ"

    dpts.units ="m"
    dpts.long_name="depth"
    dpts.standard_name ="depth"
    dpts.globwave_name="depth"
    dpts.valid_min=0.
    dpts.valid_max=10000.
    dpts.scale_factor=1.
    dpts.add_offset = 0.
    dpts.content = "TX"
    dpts.associates= "time station"
    #ds_bc['depth'].attrs['_FillValue']=9.96921e+36

    wss.units ="m s-1"
    wss.long_name="wind speed at 10m"
    wss.standard_name ="wind speed "
    wss.globwave_name="wind speed "
    wss.valid_min=0.
    wss.valid_max=100.
    wss.scale_factor=1.
    wss.add_offset = 0.
    wss.content = "TX"
    wss.associates = "time station"
    #ds_bc['u10m'].attrs['_FillValue']=9.96921e+36

    wdirs.units ="degree"
    wdirs.long_name="wind direction"
    wdirs.standard_name ="wind_from_direction"
    wdirs.globwave_name="wind_from_direction"
    wdirs.valid_min=0.
    wdirs.valid_max=360.
    wdirs.scale_factor=1.
    wdirs.add_offset= 0.
    wdirs.content = "TX"
    wdirs.associates = "time station"
    #ds_bc['udir'].attrs['_FillValue']=9.96921e+36

    curs.units ="m s-1"
    curs.long_name="sea water speed"
    curs.standard_name ="sea_water_speed"
    curs.globwave_name="sea_water_speed"
    curs.valid_min=0.
    curs.valid_max=100.
    curs.scale_factor=1.
    curs.add_offset = 0.
    curs.content = "TX"
    curs.associates = "time station"
    #ds_bc['curr'].attrs['_FillValue']=9.96921e+36

    curdirs.units ="degree"
    curdirs.long_name="direction from of sea water velocity"
    curdirs.standard_name ="direction_of_sea_water_velocity"
    curdirs.globwave_name="direction_of_sea_water_velocity"
    curdirs.valid_min=0.
    curdirs.valid_max=360.
    curdirs.scale_factor=1.
    curdirs.add_offset = 0.
    curdirs.content = "TX"
    curdirs.associates = "time station"
    #ds_bc['currdir'].attrs['_FillValue']=9.96921e+36
    return ds_bc


class LBCindex_getter():
    def __init__(self, wd, bc, parent, wind):
        self.conf = getConfigurationByID(os.path.dirname(wd), 'config')
        self.bc = bc
        self.parent = parent
        self.windFiles = wind
        self.wd = wd
        self.getIdxs()

    def flattenGrid(self, parentLon,parentLat):

        xx, yy = np.meshgrid(parentLon,parentLat)
        x = xx.flatten()
        y = yy.flatten()
        return  np.array([x, y])


    def getBCcoords(self):
        # extract all grid nodes as dict from grd file. Dict key is the id
        with open(os.path.join(self.conf.base, self.conf.grid.bottomFile), 'r') as fileIn:
            copy = False
            nodes = {}
            for line in fileIn.readlines():
                if line.strip() == "$Nodes":
                    copy = True
                elif line.strip() == "$EndNodes":
                    copy = False
                elif copy:
                    l = line.split()
                    nodes[l[0]] = l[1:]
        self.nodesDict = nodes

    def getId(self, dataset, datasetType,variableForLand, outname):
        force= self.conf.model.lateralBoundaries.refreshIndices.upper() in ["T", "TRUE"]
        if force or not os.path.exists(os.path.join(os.path.dirname(self.bc), outname)):
            ids = self.childBIDnearest(dataset[-1], datasetType,variableForLand, outname)
        else:
            ids = np.loadtxt(os.path.join(os.path.dirname(self.bc), outname))
        return ids

    def nearest1Dcoords(self, dataset, coords,datasetType, variableForLand,xy):
        """

        :param dataset: netcdf
        :param coords: x,y,z from child model
        :param datasetType: parentWave or meteoData
        :param variableForLand: variable needed to find point on sea
        :param xy: flatten coordinates from parent model
        :return:
        """
        # return array[ parentIdLon, parentIdLat, depth]
        x=float(coords[0])
        y=float(coords[1])
        depth=np.abs(float(coords[-1]))

        idxs = np.argsort(np.sqrt((xy[0] - x) ** 2 + (xy[1] - y) ** 2))
        data=dataset[variableForLand]
        if len(data.values.shape)==3:
            subSetter = {self.conf.copernicusFiles[datasetType].time: 0}

        elif len(data.values.shape)==4:
            subSetter = {self.conf.copernicusFiles[datasetType].time: 0,
                         self.conf.copernicusFiles[datasetType].depth:0}
        elif len(data.values.shape)==5:
            subSetter = {self.conf.copernicusFiles[datasetType].time: 0,
                         self.conf.copernicusFiles[datasetType].freq: 0,
                         self.conf.copernicusFiles[datasetType].dir: 0}
        data = data.isel(subSetter)
        
        xx=xy[0]
        yy=xy[1] 
        for idx in idxs:
            s={self.conf.copernicusFiles[datasetType].lon: xx[idx],
            self.conf.copernicusFiles[datasetType].lat:yy[idx]}
            near = data.sel(s, method='nearest').values
            if np.isnan(near):
                print (idx, 'skipped')
                pass
            else:
                print("Found: %s, %s" % (yy[idx], xx[idx]))
                break

        iLa,iLo=np.unravel_index(idx, shape=data.shape)
        print(data[iLa,iLo].values)
        #iLa = np.argmin(np.abs(( - dataset[latCoord].data)))
        #iLo = np.argmin(np.abs((float(coords[0]) - dataset[lonCoord].data)))
        return [iLo, iLa, depth]


    def childBIDnearest(self, dataset, datasetType,variableForLand, outname):
        """

        :param dataset: list of netcdf (wind or wave)
        :param latCoord:
        :param lonCoord:
        :param outname: is the txt filename (idAtmos or idWave)
        :return:
        """

        # use id of nodes in list of boundaries to get  child coordinates and nearest parent id
        childIds = np.loadtxt(self.bc).astype(int)
        self.getBCcoords()
        # this check is needed for sv07 waves in which we have 2 files for each day, to flat the list
        if isinstance(dataset[0], list):
            dataset = [f for ff in dataset for f in ff]
            ds = myMFdataset(dataset)
        else:
            ds = xr.open_dataset(dataset)

        xy = self.flattenGrid(ds[self.conf.copernicusFiles[datasetType].lon], ds[self.conf.copernicusFiles[datasetType].lat])
        idxs = []
        for i in childIds:
            childCoords = self.nodesDict[str(i)]
            idxs.append(self.nearest1Dcoords(ds,childCoords, datasetType,variableForLand,xy))

        out = np.array(idxs)
        np.savetxt(os.path.join(os.path.dirname(self.bc), outname), out, fmt=['%i', '%i', '%.1f'])
        return out

    def getIdxs(self):
        """

        :return: Indices array (idlatparent, idLonparent depthparent, idLatwind idlonWind)
        """
#self.conf.copernicusFiles.meteoData.lat,self.conf.copernicusFiles.meteoData.lon,self.conf.copernicusFiles.meteoData.u,
        # WIND
        windIds = self.getId(self.windFiles,'meteoData',self.conf.copernicusFiles.meteoData.u,  'idAtmos.txt')



        # WAVE
        if self.conf.model.lateralBoundaries.type =='computed':
            variable4land = self.conf.copernicusFiles.parentWave.hs
        else:
            variable4land = self.conf.copernicusFiles.parentWave.efth

        parentIds = self.getId(self.parent, 'parentWave',variable4land, 'idWave.txt')

        self.windDs = myMFdataset(self.windFiles)

        print(parentIds.shape, windIds.shape)
        uniques = duplicatesId(parentIds)

        concatId = np.hstack((parentIds[uniques], windIds[:, :2][uniques]))

        self.ids = concatId

        print(self.ids.shape)

class LBC_writer():
    def __init__(self, wd, conf, wind, parent):
        self.wd = wd
        self.wind = wind
        self.parent = parent
        self.conf = getConfigurationByID(conf, 'config')
        self.spectraDefinition()

    def spectraDefinition(self):
        """

        :return: add to the LBCwriter class the frequencies and the directions in RAD
        """
        self.freqs = getFreqs(self.conf.grid.minFrequency, self.conf.grid.NoFrequencies,
                              self.conf.grid.frequencyIncrement)

        self.dirs_rad = degToRad(getDirs(self.conf.grid.NoDirections))

        self.dirs_deg = getDirs(self.conf.grid.NoDirections)

    def bc_header(self):
        """

        :return: write the header for each boundary TEXT file
        """
        outname = os.path.join(self.wd, 'bc_header.spc')

        fbins = [np.format_float_scientific(f, precision=3, exp_digits=2, unique=False) for f in self.freqs]
        f_chunks = (' '.join(str(l) + '\n' * (n % 8 == 7) for n, l in enumerate(fbins)))

        dbins = [np.format_float_scientific(f, precision=3, exp_digits=2, unique=False) for f in self.dirs_rad]
        d_chunks = (' '.join(str(l) + '\n' * (n % 8 == 7) for n, l in enumerate(dbins)))

        if d_chunks[-1] != '\n':
            d_chunks = d_chunks + '\n'

            with open(outname, 'w')as o:
                o.write("'WAVeWATCH III SPeCTRA'     %02d    %02d     1 'MULTIGRID GLOBAL05+OTHeRS     '\n " % (
                    self.conf.grid.NoFrequencies,
                    self.conf.grid.NoDirections))
                o.write(f_chunks)
                o.write('\n ')
                o.write(d_chunks)
                o.close()

    def getWindData(self, point):
        subSetter = {self.conf.copernicusFiles.meteoData.lon: int(point[0]),
                     self.conf.copernicusFiles.meteoData.lat: int(point[1])}
        windData = self.windDs.isel(subSetter)
        vVar = windData[self.conf.copernicusFiles.meteoData.v]
        uVar = windData[self.conf.copernicusFiles.meteoData.u]

        wDir = (np.arctan2(vVar, uVar) * 180. / np.pi)
        wSpeed = np.sqrt((uVar ** 2.) + (vVar ** 2.))

        return wDir, wSpeed

    def onePointBC_netcdf(self, ids, pointName):
        if self.conf.model.lateralBoundaries.type =='computed':
            computed = True
        else:
            computed = False
        spName = 'id_{}'.format(pointName)
        outFile = os.path.join(self.wd, '{}_spec.nc'.format(spName))

        waveLatId = int(ids[1])
        waveLonId = int(ids[0])
        dpt = float(ids[2])
        wndId = list(map(int, [ids[3], ids[4]]))

        # WAVE PARAMETERS
        subSetter = {self.conf.copernicusFiles.parentWave.lat: waveLatId,
                     self.conf.copernicusFiles.parentWave.lon: waveLonId
                     }

        self.windDs = myMFdataset(self.wind)
        parMod = myMFdataset(self.parent).isel(subSetter)

        hss = parMod[self.conf.copernicusFiles.parentWave.hs].load()
        tps = parMod[self.conf.copernicusFiles.parentWave.tp].load()
        dirs = parMod[self.conf.copernicusFiles.parentWave.dir].load()
        # if self.conf.model.lateralBoundaries.computed.forced.flag.upper() in ["T", "TRUE"]:
        #     hss.values=np.zeros_like(hss)+ float(self.conf.model.lateralBoundaries.computed.forced.hs)
        #     tps.values=np.zeros_like(tps)+float(self.conf.model.lateralBoundaries.computed.forced.tp)
        #     dirs.values=np.zeros_like(dirs)+float(self.conf.model.lateralBoundaries.computed.forced.dir)
        # else:
        for i, t in enumerate(hss):
            hss[i] = seaoverland(np.ma.masked_array(t, mask=np.isnan(t)), 10)

        for i, t in enumerate(tps):
            tps[i] = seaoverland(np.ma.masked_array(t, mask=np.isnan(t)), 10)

        for i, t in enumerate(dirs):
            dirs[i] = seaoverland(np.ma.masked_array(t, mask=np.isnan(t)), 10)

        dirs = dirs.values
        nDir = self.conf.grid.NoDirections
        lo = np.round(parMod[self.conf.copernicusFiles.parentWave.lon].values, 2)
        la = np.round(parMod[self.conf.copernicusFiles.parentWave.lat].values, 2)

        # WIND PARAMETERS
        windDirs, windSpeeds = self.getWindData(wndId)

        ds_bc=getBCdataset(parMod.time, pointName, self.freqs, self.dirs_deg, outFile, computed=computed)
        buffer_input = []  # this is just to save what is used as input to rebuild spectra
        for i,t in enumerate(parMod.time):
            print (i)
            wSpeed = windSpeeds.sel({self.conf.copernicusFiles.meteoData.time: t}, method='nearest').values
            wDir = windDirs.sel({self.conf.copernicusFiles.meteoData.time: t}, method='nearest').values

            hs = hss.sel({self.conf.copernicusFiles.parentWave.time: t}).values
            tp = tps.sel({self.conf.copernicusFiles.parentWave.time: t}).values
            input_dir = dirs[i]
            conv_dir = fixBCdir(input_dir)
            specifics = f"time_{t.values},coords_y{la} x{lo}, hs_{hs},tp_{tp}, dir_{input_dir}_{conv_dir}, wsp_{wSpeed}, wdir_{wDir}"
            buffer_input.append(specifics)
            nonDimSpec, dimSpec = JONSWAP(hs, tp, conv_dir, nDir, self.freqs).main()
            print (dimSpec.T.shape, ds_bc['efth'].shape)
            ds_bc['longitude'][i]=  lo
            ds_bc['latitude'][i] = la
            ds_bc['depth'][i] = dpt
            ds_bc['u10m'][i]=  wSpeed
            ds_bc['udir'][i] = wDir
            ds_bc['efth'][i] = dimSpec.T

        #with open(os.path.join(self.wd, '{}_specs.txt'.format(spName)), 'w') as f:
        #    f.write('\n'.join(buffer_input))

        ds_bc.close()


    def onePointBC(self, ids, pointName):
        self.bc_header()

        waveLatId = int(ids[1])
        waveLonId = int(ids[0])
        dpt = float(ids[2])

        wndId = list(map(int, [ids[3], ids[4]]))
        spName = 'id_{}'.format(pointName)
        outFile = os.path.join(self.wd, '{}.spc'.format(spName))
        # WAVE PARAMETERS
        subSetter = {self.conf.copernicusFiles.parentWave.lat: waveLatId,
                     self.conf.copernicusFiles.parentWave.lon: waveLonId
                     }
        self.windDs = xr.open_mfdataset(self.wind, combine='by_coords')
        parMod = xr.open_mfdataset(self.parent, combine='by_coords')

        hss = parMod[self.conf.copernicusFiles.parentWave.hs].load()
        tps = parMod[self.conf.copernicusFiles.parentWave.tp].load()
        dirs = parMod[self.conf.copernicusFiles.parentWave.dir].load()

        for i, t in enumerate(hss):
            hss[i] = seaoverland(np.ma.masked_array(t, mask=np.isnan(t)), 10)
        hss = hss.isel(subSetter)
        for i, t in enumerate(tps):
            tps[i] = seaoverland(np.ma.masked_array(t, mask=np.isnan(t)), 10)
        tps = tps.isel(subSetter)
        for i, t in enumerate(dirs):
            dirs[i] = seaoverland(np.ma.masked_array(t, mask=np.isnan(t)), 10)
        dirs = dirs.isel(subSetter)

        nDir = self.conf.grid.NoDirections
        lo = np.round(parMod[self.conf.copernicusFiles.parentWave.lon].data[waveLonId], 2)
        la = np.round(parMod[self.conf.copernicusFiles.parentWave.lat].data[waveLatId], 2)

        # WIND PARAMETERS
        wDirs, wSpeeds = self.getWindData(wndId)

        with open(outFile, 'w') as o:

            with open(os.path.join(self.wd, 'bc_header.spc')) as head:
                for line in head:
                    o.write(line)
                head.close()

            for t in parMod.time:

                wSpeed = wSpeeds.sel({self.conf.copernicusFiles.meteoData.time: t}, method='nearest').data
                wDir = wDirs.sel({self.conf.copernicusFiles.meteoData.time: t}, method='nearest').data
                timestamp = datetime.strptime(str(t.data).split('.')[0], "%Y-%m-%dT%H:%M:%S").strftime(
                    '%Y%m%d %H%M%S')

                hs = hss.sel({self.conf.copernicusFiles.parentWave.time: t})
                tp = tps.sel({self.conf.copernicusFiles.parentWave.time: t})
                dir = dirs.sel({self.conf.copernicusFiles.parentWave.time: t})

                nonDimSpec, dimSpec = JONSWAP(hs, tp, dir, nDir, self.freqs).main()
                print('hs:', hs.data, ' dir:', dir.data, ' coords:', lo, la, 'time:', t.data)
                o.write('%s\n' % timestamp)
                o.write("'{pointId}'{lat}{lon}{depth}{wSpeed}{wDir}{cSpeed}{cDir}\n  ".format(
                    pointId=stringAllocation(10, spName[:10], 'forth'),
                    lat='%7.2f' % la,
                    lon='%7.2f' % lo,
                    depth='%10.1f' % float(dpt),
                    wSpeed='%7.2f' % np.round(wSpeed, 2),
                    wDir='%6.1f' % np.round(wDir, 1),
                    cSpeed='%7.2f' % 0.0,
                    cDir='%6.1f' % 270.0))

                spc = [np.format_float_scientific(i, precision=3, exp_digits=2, unique=False) for i in
                       dimSpec.ravel()]
                chunkSpc = ('  '.join(str(l) + '\n' * (n % 7 == 6) for n, l in enumerate(spc)))
                if chunkSpc[-1] != '\n':
                    chunkSpc = chunkSpc + '\n'
                o.write(chunkSpc)
            o.write('\n')
        o.close()


class LBC_extractor():
    def __init__(self, wd, conf, wind, parent):
        self.wd = wd
        self.windDs = wind
        self.parMod = parent
        self.conf = getConfigurationByID(conf, 'config')
        self.spectraDefinition()


    def spectraDefinition(self):
        """

        :return: add to the LBCwriter class the frequencies and the directions in RAD
        """

        if self.conf.model.lateralBoundaries.type.upper() in ['FROMPARENT']:
            self.freqs=self.parMod[self.conf.copernicusFiles.parentWave.freq].values
            self.dirs_deg=self.parMod[self.conf.copernicusFiles.parentWave.dir].values

        else:
            self.freqs = getFreqs(self.conf.grid.minFrequency, self.conf.grid.NoFrequencies,
                                  self.conf.grid.frequencyIncrement)

            self.dirs_rad = degToRad(getDirs(self.conf.grid.NoDirections))

            self.dirs_deg = getDirs(self.conf.grid.NoDirections)



    def getWindData(self, point):
        subSetter = {self.conf.copernicusFiles.meteoData.lon: int(point[0]),
                     self.conf.copernicusFiles.meteoData.lat: int(point[1])}
        windData = self.windDs.isel(subSetter)
        vVar = windData[self.conf.copernicusFiles.meteoData.v]
        uVar = windData[self.conf.copernicusFiles.meteoData.u]

        wDir = (np.arctan2(vVar, uVar) * 180. / np.pi)  # +90
        wSpeed = np.sqrt((uVar ** 2.) + (vVar ** 2.))

        return wDir, wSpeed

    def getSpectraBC(self, ids, pointName):

        spName = 'id_{}'.format(pointName)
        outFile = os.path.join(self.wd, '{}_spec.nc'.format(spName))

        waveLonId = int(ids[0])
        waveLatId = int(ids[1])
        dpt = float(ids[2])
        wndId = list(map(int, [ids[3], ids[4]]))

        # WAVE PARAMETERS
        subSetter = {self.conf.copernicusFiles.parentWave.lat: waveLatId,
                     self.conf.copernicusFiles.parentWave.lon: waveLonId
                     }
        efths = self.parMod[self.conf.copernicusFiles.parentWave.efth]
        efths = efths.isel(subSetter)

        lo = self.parMod[self.conf.copernicusFiles.parentWave.lon].data[waveLonId]
        la = self.parMod[self.conf.copernicusFiles.parentWave.lat].data[waveLatId]

        # WIND PARAMETERS
        windDirs, windSpeeds = self.getWindData(wndId)

        ds_bc=getBCdataset(self.parMod.time, pointName, self.freqs, self.dirs_deg, outFile)

        for i,t in enumerate(self.parMod.time):
            print (i)
            wSpeed = windSpeeds.sel({self.conf.copernicusFiles.meteoData.time: t}, method='nearest').values
            wDir = windDirs.sel({self.conf.copernicusFiles.meteoData.time: t}, method='nearest').values
            efth = efths.sel({self.conf.copernicusFiles.parentWave.time: t}).values
            #print (efth)
            ds_bc['longitude'][i]=  lo
            ds_bc['latitude'][i] = la
            ds_bc['depth'][i] = dpt
            ds_bc['u10m'][i]=  wSpeed
            ds_bc['udir'][i] = wDir
            ds_bc['efth'][i] = efth
            print (ds_bc['efth'][i][0].shape)
        ds_bc.close()


