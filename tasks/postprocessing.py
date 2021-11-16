import shutil
from glob import glob
from netCDF4 import Dataset
import numpy as np
import os
import xarray as xr
from glob import glob
import subprocess
from natsort import natsorted
from utils import checkFunctionCompleted
import datetime



class Postprocessing:
    def __init__(self, conf, paths,args,checks,submitter):
        self.conf=conf
        self.checks=checks
        self.args=args
        self.dryrun=args.dryRun
        if str(self.dryrun).upper() in ['Y','YES']:
            self.dryRun=True
        else:
            self.dryRun = False
        if not self.checks.NCcomplete():
            self.getNetCDF(submitter)
        self.restoreRundir(paths)

    def getNetCDF(self,submitter):
        if not self.dryRun:
            submitter.systemCommand("{}".format(os.path.join(self.conf.model.executable, 'ww3_ounf')))
            open('netCDF_conversionComplete', 'w').close()

    def restoreRundir(self, paths):
        if not self.dryRun:
            ncList=natsorted(glob(os.path.join(paths.runDir, 'ww3.*.nc')))

            #OPERATIONAL RUN
            # in op run ncList =ncList[1:]
            for nc in ncList:
                fname=os.path.basename(nc)
                shutil.move(nc,os.path.join(paths.ncDir,fname))
            dayRst=os.path.join(paths.workingDir,'restart_{}.ww3'.format(self.args.startDate))
            newNameRst=os.path.join(paths.runDir,'restart_{}.ww3'.format(self.args.startDate))
            lastRst=os.path.join(paths.runDir, 'restart001.ww3')
            print ('last restart available: ',lastRst)
            print ('renamed: to',newNameRst)
            print ('copied to ', dayRst)
            shutil.move(lastRst, newNameRst)
            shutil.copy(newNameRst, dayRst)

        os.chdir(os.path.join(self.conf.base, 'tasks'))
        shutil.move(paths.runDir, os.path.join(paths.runsArchive,paths.runDir.split('/')[-1]))


class Regridding:
    def __init__(self,conf,paths):
        self.conf=conf
        self.paths=paths

    def regriddingMain(self, area):
        os.chdir(self.paths.ncDir)
        self.regridOutput(area)
        self.mergeRegrid( area)

    def regridOutput(self, area):
        print('regridding output')
        vars = self.conf.post.regrid.area[area].variables
        res = self.conf.post.regrid.area[area].resolution

        mask = os.path.join(os.path.dirname(self.conf.grid.bottomFile),  'mask%s_%s.nc' % (int(res * 100000), area))
        ncs = natsorted([os.path.join(self.paths.ncDir, f) for f in os.listdir(self.paths.ncDir) if f.endswith('.nc') if
               f.startswith('ww3.')])

        box = (',').join([str(f) for f in self.conf.post.regrid.area[area].boundingBox]).replace('[', '').replace(']','')

        if not os.path.exists(mask):
            subprocess.call(
                'python {exe} rast mask {nc} {msk} --dx {res} --dy {res} --topology {topo} --box {bbox}'.format(
                    bbox=box, exe=self.conf.post.regrid.regridderPath,
                    nc=ncs[0], msk=mask, res=res,
                    topo=self.conf.post.regrid.topology), shell=True)
        n = 0
        startTs = 0
        for nc in ncs:
            n += len(xr.open_dataset(nc).time)

            if area == 'domain':
                bbox = box.replace(',', '_')
            else:
                bbox = area

            for var in vars:
                outfile = 'REG_{startTs}_{n}_0_{var}_{bbox}_{res}.nc_reg'.format(startTs=startTs, n=n, var=var,
                                                                                 bbox=bbox, res=res)

                cmd = 'python {exe} rast var {ncs} {msk}  {out} --var {var}'.format(
                    exe=self.conf.post.regrid.regridderPath,
                    ncs=nc, msk=mask, out=outfile, var=var)

                subprocess.call(cmd, shell=True)
            startTs = n

    def mergeRegrid(self, area):
        print('merging regridded files')
        box = (',').join([str(f) for f in self.conf.post.regrid.area[area].boundingBox]).replace('[', '').replace(']',
                                                                                                                  '')
        vars = self.conf.post.regrid.area[area].variables
        if area == 'domain':
            bbox = box.replace(',', '_')
        else:
            bbox = area
        for var in vars:
            ncs = natsorted(
                [os.path.join(self.paths.ncDir, f) for f in os.listdir(self.paths.ncDir) if f.endswith('.nc_reg') if
                 f.__contains__(var) if f.__contains__(bbox)])
            print(ncs)
            outfile = 'REG_{bbox}_{var}.nc'.format(var=var, bbox=bbox)

            if not os.path.exists(os.path.join(self.paths.ncDir, outfile)):
                out = xr.open_mfdataset(ncs)
                out.to_netcdf(outfile)
            else:
                pass
