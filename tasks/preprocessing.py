import os
import shutil
import subprocess
from glob import glob
from natsort import natsorted
from boundaryConditions import LBCindex_getter,LBC_writer,LBC_extractor
from inputFieldProvider import buildFilePath, findLastForecast, forecastAnalysisPath, WNSgetter, WNDgetter, \
    CURgetter, LEVgetter,selectInputFile
from templates import Templates
from utils import checkFunctionCompleted,getDaysBetweenDates, myMFdataset
from bjobs import Bjobs
import numpy as np
import xarray as xr


class Preprocessing():
    def __init__(self, conf,rundir,runDuration,args, checks,submitter,lastrun):
        self.checks = checks
        self.conf = conf
        self.rundir = rundir
        self.startdate = args.startDate
        self.workingDir = args.workingDir
        self.submitter=submitter
        self.runDuration=runDuration
        self.lastrun=lastrun

    def computeLBC(self):
        print ('LBC processing...')
        LBC = LateralBoundaryCondition(self.conf, self.rundir,self.workingDir, self.startdate,self.runDuration).compute()
        open('LBC_complete', 'w').close()

    def extractLBC(self):
        print ('LBC processing...')
        LBC = LateralBoundaryCondition(self.conf, self.rundir,self.workingDir, self.startdate,self.runDuration).extract()
        open('LBC_complete', 'w').close()

    def computeSBC(self):
        print ('SBC processing...')
        SBC = SurfaceBoundaryCondition(self.conf, self.rundir, self.startdate, self.runDuration,self.workingDir)
        SBC.main()
        for field in SBC.fields:
            # iterate prep executable on all fields
            # PREP
            #SBC.writePrepINP(field)
            #self.submitter.systemCommand("{}".format(os.path.join(self.conf.model.executable, 'ww3_prep')))
            #PRNC
            SBC.writePrncINP(field)
            self.submitter.systemCommand("{}".format(os.path.join(self.conf.model.executable, 'ww3_prnc')))
        open('SBC_complete', 'w').close()

    def processTemplates(self):
        #try:
        if self.conf.model.lateralBoundaries.type.upper() in ['FROMPARENT', 'COMPUTED'] :
            #spectra= natsorted(glob(os.path.join(self.rundir, 'id_*.spc')))
            spectra = natsorted(glob(os.path.join(self.rundir, 'id_*spec.nc')))
            spectra_def=xr.open_dataset(spectra[0])
        else:
            spectra=False
            spectra_def=False
        Templates(self.conf, self.startdate,self.runDuration, self.rundir,spectra,spectra_def,self.lastrun).generate()
        #    open('TMPL_complete', 'w').close()
        #except:
        #    print ('error in Templates')

    def processINP(self):
        outgrid=os.path.join(self.workingDir,'mod_def.ww3')
        if not os.path.exists(outgrid):
            self.submitter.systemCommand("{}".format(os.path.join(self.conf.model.executable, 'ww3_grid')))
            shutil.copy(os.path.join(self.rundir,'mod_def.ww3'),outgrid)
        else:
            shutil.copy(outgrid, os.path.join(self.rundir,'mod_def.ww3'))



    def main(self):

        print ('preprocessing starting')
        if self.conf.model.lateralBoundaries.type.upper() not in ['COMPUTED','FROMPARENT','F', 'FALSE']:
            exit('please check input at model.lateralBoundaries.type in conf.yaml. Valid keys are: computed, fromParent, F')
        os.chdir(self.rundir)
        if self.conf.model.lateralBoundaries.type =='computed':
            print ('computing boundary conditions')
            checkFunctionCompleted(self.checks.LBCcomplete,self.computeLBC)
        if self.conf.model.lateralBoundaries.type =='fromParent':
            print ('extracting boundary conditions')
            #exit('LBC from parent model not implemented yet')
            checkFunctionCompleted(self.checks.LBCcomplete,self.extractLBC)
        print ('filling templates')
        checkFunctionCompleted(self.checks.TMPLcomplete,self.processTemplates)
        print ('getting symlink to grid')
        self.getSimlink()
        print ('preparing .inp files')
        checkFunctionCompleted(self.checks.INP_gridComplete, self.processINP)
        print ('computing surface boundary condition')
        checkFunctionCompleted(self.checks.SBCcomplete,self.computeSBC)
        if self.conf.model.lateralBoundaries.type.upper() in ['COMPUTED', 'FROMPARENT']:
            self.submitter.systemCommand(os.path.join(self.conf.model.executable, 'ww3_bounc'))
            #self.submitter.systemCommand(os.path.join(self.conf.model.executable, 'ww3_bound')) 
        print ('restart managing')
        Restart(self.conf,self.workingDir,self.rundir).manage(self.submitter)

    def getSimlink(self):
        # Create grid files symlink to run path
        if self.conf.grid.type == 'UNST':
            simlinks = [os.path.join(self.conf.base, self.conf.grid.bottomFile),
                        os.path.join(self.conf.base, self.conf.grid.obstFile_local),
                        os.path.join(self.conf.base, self.conf.grid.obstFile_shadow)]
            #TSG
            #simlinks = [os.path.join(self.conf.base, self.conf.grid.bottomFile),os.path.join(self.conf.base,'data/msh/sgmask.txt')]
             
        else:
            simlinks = [os.path.join(self.conf.base, self.conf.grid.bottomFile),
                        os.path.join(self.conf.base, self.conf.grid.obstFile_regular),
                        os.path.join(self.conf.base, self.conf.grid.maskFile)]
        for f in simlinks:
            try:
                os.symlink(f, os.path.join(self.rundir, os.path.basename(f)))
            except:
                pass


class LateralBoundaryCondition():
    def __init__(self,conf,rundir, workingdir,startdate,duration):
        self.conf=conf
        self.rundir=rundir
        self.workingdir=workingdir
        self.startdate=startdate
        self.runDuration=duration


    def compute(self):
        print('Computing Indexes for LBC')
        days=getDaysBetweenDates(self.startdate,self.runDuration+1)

        windFiles = []
        for day in days:
            wnds = selectInputFile(self.conf, 'meteoData', day)
            if isinstance(wnds, list):
                for wnd in wnds:
                    if wnd not in windFiles:
                        windFiles.append(wnd)
            else:
                windFiles.append(wnds)

        #windFile=[selectInputFile(self.conf,'meteoData',day)for day in days]
        parentFile = []
        for day in days:
            prnts = selectInputFile(self.conf, 'parentWave', day)
            if isinstance(prnts, list):
                for prnt in prnts:
                    if prnt not in parentFile:
                        parentFile.append(prnt)
            else:
                parentFile.append(prnts)

        #parentFile = [selectInputFile(self.conf, 'parentWave', day) for day in days]

        # wd, conf, bc, parent, wind
        lbc = LBCindex_getter(self.rundir,
                              os.path.join(self.conf.base,self.conf.model.lateralBoundaries.path),
                              parentFile,
                              windFiles)
        allCMD=[]
        for i, id in enumerate(lbc.ids):
            i = i + 1
            strId = (' ').join([str(j) for j in id.tolist()])
            argv = [strId, i, self.rundir, self.workingdir, days]
            print(argv)
            #allCMD.append('bsub -P {p} -q {q} -J {i}_bc python {executable} {argv}'.format(
            #        executable=os.path.join(self.conf.base, 'tasks/runOneLBC.py'),p=self.conf.model.project_queue, i=i, q=self.conf.model.serial_queue, argv=argv))
            allCMD.append('bsub  -P {p} -q {q} -M {m} -R rusage[mem={m}] -J {i}_bc python {executable} {argv}'.format(
                    executable=os.path.join(self.conf.base, 'tasks/runOneLBC_netcdf.py'),p=self.conf.model.project_queue, i=i,m=self.conf.model.lateralBoundaries.memoryUsage, q=self.conf.model.serial_queue, argv=argv))
            if i == len(lbc.ids):
                Bjobs(allCMD,self.conf.model.project_queue)
                allCMD = []
            else:
                pass

            #if i == len(lbc.ids):
            #    open(os.path.join(self.rundir, 'bc.done'), 'w').close()

        self.spectra= natsorted(glob(os.path.join(self.rundir, 'id_*.spc')))

    def getIDXSboxes(self,idxs):
        """

        :param idxs: indices from LBCindex_getter
        :return: parent, wind arrays with box indices to cut parent and wind files. Each of them is [minX,maxX,minY,maxY]. The size of the box is +2 for X and Y
        """

        latParent=idxs[:,1]
        lonParent = idxs[:, 0]
        latWind = idxs[:, 4]
        lonWind = idxs[:,3]

        boxParent=[np.nanmin(lonParent)-2,np.nanmax(lonParent)+2,np.nanmin(latParent)-2,np.nanmax(latParent)+2]
        boxWind = [np.nanmin(lonWind)-2, np.nanmax(lonWind)+2, np.nanmin(latWind)-2, np.nanmax(latWind)+2]

        return boxParent,boxWind

    def converterBoxIdxs(self,ids):
        ids.T[0]-=np.nanmin(ids.T[0])
        ids.T[1]-=np.nanmin(ids.T[1])
        ids.T[3]-=np.nanmin(ids.T[3])
        ids.T[4]-=np.nanmin(ids.T[4])
        
        ids.T[0]+=2
        ids.T[1]+=2
        ids.T[3]+=2
        ids.T[4]+=2
        return ids

    def extract(self):
        print('Extracting Indexes for LBC')
        days = getDaysBetweenDates(self.startdate, self.runDuration + 1)

        windFiles = []
        for day in days:
            wnds = selectInputFile(self.conf, 'meteoData', day)
            if isinstance(wnds, list):
                for wnd in wnds:
                    if wnd not in windFiles:
                        windFiles.append(wnd)
            else:
                windFiles.append(wnds)

        #windFiles = [selectInputFile(self.conf, 'meteoData', day) for day in days]

        spectraFiles = []
        for day in days:
            spcs = selectInputFile(self.conf, 'parentWave', day)
            if isinstance(spcs, list):
                for spc in spcs:
                    if spc not in spectraFiles:
                        spectraFiles.append(spc)
            else:
                spectraFiles.append(spcs)

        #spectraFiles = [selectInputFile(self.conf, 'parentWave', day) for day in days]


        # wd, conf, bc, parent, wind
        lbc = LBCindex_getter(self.rundir,
                              os.path.join(self.conf.base, self.conf.model.lateralBoundaries.path),
                              spectraFiles,
                              windFiles)

        boxParent, boxWind=self.getIDXSboxes( lbc.ids)

        winds=myMFdataset(windFiles,boxWind)
        spectra=myMFdataset(spectraFiles,boxParent)


        writer=LBC_extractor(self.rundir, self.workingdir, winds, spectra)

        for i, id in enumerate(self.converterBoxIdxs(lbc.ids)):
            i = i + 1
            strId = (' ').join([str(j) for j in id.tolist()])
            argv = [strId, i,  days]
            print(argv)
            # allCMD.append('bsub -P {p} -q {q} -J {i}_bc python {executable} {argv}'.format(
            #        executable=os.path.join(self.conf.base, 'tasks/runOneLBC.py'),p=self.conf.model.project_queue, i=i, q=self.conf.model.serial_queue, argv=argv))
            #allCMD.append('bsub -P {p} -q {q} -J {i}_bc python {executable} {argv}'.format(
            #    executable=os.path.join(self.conf.base, 'tasks/runOneLBC_netcdf.py'),
            #    p=self.conf.model.project_queue, i=i, q=self.conf.model.serial_queue, argv=argv))
            writer.getSpectraBC(id, i)

        self.spectra = natsorted(glob(os.path.join(self.rundir, 'id_*.spc')))




class SurfaceBoundaryCondition():
    def __init__(self,conf, rundir, startdate,runDuration,workingDir):
        self.conf=conf
        self.rundir=rundir
        self.startdate=startdate
        self.workingDir=workingDir
        self.runDuration = runDuration

        #self.main()

    def main(self):

        print('starting with preprocessing Surface Boundaries')
        self.fields=[]
        if str(self.conf.copernicusFiles.meteoData.ww3Name).upper() == 'WNS':
            WNSgetter(self.conf, self.startdate, self.runDuration, self.conf.grid.fillValue).getWNS()
            self.fields.append('WNS')
            print('WNS processed')
        if (str(self.conf.copernicusFiles.meteoData.ww3Name).upper() == 'WND')&(str(self.conf.model.forcings.windFlag).upper() == 'T'):
            WNDgetter(self.conf, self.startdate, self.runDuration, self.conf.grid.fillValue).getWND()
            self.fields.append('WND')
            print('WND processed')
        else:
            print('No WIND')
        if str(self.conf.model.forcings.curFlag).upper() == 'T':
            CURgetter(self.conf, self.startdate, self.runDuration, self.conf.grid.fillValue).getCUR()
            self.fields.append('CUR')
            print('CUR processed')
        if str(self.conf.model.forcings.levFlag).upper() == 'T':
            LEVgetter(self.conf, self.startdate, self.runDuration, self.conf.grid.fillValue).getLEV()
            self.fields.append('LEV')
            print('LEV processed')

    def writePrepINP(self, field):
        with open(os.path.join(self.rundir, 'ww3_prep.inp'), 'w')as out:
            spec = os.path.join(self.rundir, '{}_gridSpecs.txt'.format(field))
            with open(spec, 'r') as s:
                out.write('$\n ')
                for line in s:
                    out.write(line)
                out.write('\n')
                out.write("'NAME' 2 1 '(..T..)' '(..F..)'"
                          "\n 20 '{}.data'".format(field))

    def writePrncINP(self, field):
        with open(os.path.join(self.rundir, 'ww3_prnc.inp'), 'w')as out:
            out.write('$\n ')
            out.write(f"'{field}' 'LL' T T\n")
            out.write('$\n ')

            if field=='WND':
                out.write("{lon} {lat}\n".format(lon=self.conf.copernicusFiles.meteoData.lon,
                                                 lat=self.conf.copernicusFiles.meteoData.lat))
                out.write('$\n ')
                out.write("{U} {V}\n".format(U=self.conf.copernicusFiles.meteoData.u,V=self.conf.copernicusFiles.meteoData.v))
            elif field=='WNS':
                out.write("{lon} {lat}\n".format(lon=self.conf.copernicusFiles.meteoData.lon,
                                                 lat=self.conf.copernicusFiles.meteoData.lat))
                out.write('$\n ')
                out.write("{U} {V} DT\n".format(U=self.conf.copernicusFiles.meteoData.u,V=self.conf.copernicusFiles.meteoData.v))
            elif field=='CUR':
                out.write("{lon} {lat}\n".format(lon=self.conf.copernicusFiles.parentHydro.lon,
                                                 lat=self.conf.copernicusFiles.parentHydro.lat))
                out.write('$\n ')
                out.write("{U} {V} \n".format(U=self.conf.copernicusFiles.parentHydro.waterVelocity.u,V=self.conf.copernicusFiles.parentHydro.waterVelocity.v))
            elif field=='LEV':
                out.write("{lon} {lat}\n".format(lon=self.conf.copernicusFiles.parentHydro.lon,
                                                 lat=self.conf.copernicusFiles.parentHydro.lat))
                out.write('$\n ')
                out.write("{SSH} \n".format(SSH=self.conf.copernicusFiles.parentHydro.waterLevel.ssh))
            out.write('$\n ')
            out.write(f"'{field}.nc'")


class Restart():
    def __init__(self, conf,workingdir,rundir):
        self.workingdir=workingdir
        self.rundir=rundir
        self.conf = conf

    def manage(self,submitter):
        if self.checkPreviousRST(self.workingdir,self.rundir):
            pass
        else:
            submitter.systemCommand("{}".format(os.path.join(self.conf.model.executable,'ww3_strt')))

    def checkPreviousRST(self,workingdir,rundir):
        previousRSTs = natsorted(glob(os.path.join(workingdir, 'restart_*.ww3')))#[-1]
        if len(previousRSTs)>0:
            if rundir.split('/')[-2] != rundir.split('/')[-1]:
                print ('Hot restart')
                shutil.move(previousRSTs[-1], os.path.join(rundir, 'restart.ww3'))
                return True
            else:
                print ('Cold restart')
                return False
        else:
                print ('Cold restart')
                return False


