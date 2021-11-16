import os
from datetime import timedelta
import pybars
from utils import dataFormatter
import numpy as np


class Templates():
    def __init__(self,conf,startdate,duration,rundir,spectra,spectra_def,lastrun):
        self.conf=conf
        self.startdate=startdate
        self.rundir=rundir
        self.spectra=spectra
        self.compiler = pybars.Compiler()
        self.runDuration=duration
        self.lastrun=lastrun
        self.data=self.getDict(conf,spectra_def)

    def getDict(self,conf,spectra_def):

        startRunDate = dataFormatter(self.startdate).strftime('%Y%m%d %H%M%S')

        endRunDate = (dataFormatter(self.startdate) + timedelta(days=self.runDuration) + timedelta(
            seconds=conf.model.output.timeInterval)).strftime('%Y%m%d %H%M%S')
        outTSs=86400 / conf.model.output.timeInterval * self.runDuration

        if self.conf.model.lateralBoundaries.type.upper() in ['COMPUTED', 'F', "FALSE"]:
            frequencyIncrement=conf.grid.frequencyIncrement
            minFrequency=conf.grid.minFrequency
            NoFrequencies=conf.grid.NoFrequencies
            NoDirections=conf.grid.NoDirections
        else:
            freqs=spectra_def[conf.copernicusFiles.parentWave.freq].values
            frequencyIncrement=np.round(freqs[1]/freqs[0],2)
            minFrequency=np.min(freqs)
            NoFrequencies=len(freqs)
            NoDirections=len(spectra_def[conf.copernicusFiles.parentWave.dir].values)

        data = {'model': conf.modelName,
                'namelist': [{'n': i} for i in conf.grid.namelist],
                'outvars': conf.model.output.fields,
                'frequencyIncrement': frequencyIncrement,
                'minFrequency': minFrequency,
                'NoFrequencies': NoFrequencies,
                'NoDirections': NoDirections,
                'gridType': conf.grid.type,
                'depthLimit': conf.grid.depthLimit,
                'minDepth': conf.grid.minDepth,
                'itype': conf.model.forcings.itype,
                'inputField': [{'field': conf.copernicusFiles.parentHydro.waterVelocity.ww3Name},
                               {'field': conf.copernicusFiles.meteoData.ww3Name},
                               {'field': conf.copernicusFiles.parentHydro.waterLevel.ww3Name}],
                'global_TS': conf.timeStep.global_TS,
                'CFL_TS': conf.timeStep.CFL_TS,
                'refraction_TS': conf.timeStep.refraction_TS,
                'minimum_TS': conf.timeStep.minimum_TS,
                'timeInterval': conf.model.output.timeInterval,
                'startDate': startRunDate,
                'endDate': endRunDate,
                'outTSs': outTSs+1 if self.lastrun else outTSs,
                'bottomFile': os.path.basename(conf.grid.bottomFile),
                'curFlag': conf.model.forcings.curFlag,
                'levFlag': conf.model.forcings.levFlag,
                'windFlag': conf.model.forcings.windFlag,
                'bc_file': os.path.join(conf.base, conf.model.lateralBoundaries.path),
                'scaleFactor': conf.grid.depthScaleFactor,
                'restart': 86400* self.runDuration,
                'cores': conf.model.multiprocess.cores,
                'queue': conf.model.multiprocess.queue,
                'shell': os.path.join(conf.model.executable, 'ww3_shel'),
                'wd': self.rundir,
                'project': conf.model.project_queue,
                #'sla': conf.model.service_class
                }

        if conf.grid.type == 'UNST':
            self.ext= '.unst'
            data['noPoints'] = conf.grid.noNodes
            data['obstFile_local'] = os.path.basename(conf.grid.obstFile_local)
            data['obstFile_shadow'] = os.path.basename(conf.grid.obstFile_shadow)

        else:
            self.ext = '.reg'
            data['lonResolution'] = conf.grid.lonResolution
            data['latResolution'] = conf.grid.latResolution
            data['minLon'] = conf.grid.minLon
            data['minLat'] = conf.grid.minLat
            data['lonSize'] = conf.grid.lonSize
            data['latSize'] = conf.grid.latSize
            data['maskFile'] = os.path.basename(conf.grid.maskFile)
            data['obstFile_regular'] = os.path.basename(conf.grid.obstFile_regular)

        boundNodes = []
        try:
            fo = open(data['bc_file']).readlines()
            boundNodes = [{'nodes': '  %s 1 F' % str(i).replace('\n', '')} for i in fo]
        except Exception as e:
            print(e)

        data['bc'] = boundNodes
        data['spcs'] = []
        try:
            [data['spcs'].append({'spc': spc}) for spc in self.spectra]
        except:
            pass
        return data

    def generate(self):
        listTemplates = [f for f in os.listdir(self.conf.model.templates) if f.endswith(self.ext) or f.endswith('.tmpl')]
        [self.fillTemplate( tmpl)for tmpl in listTemplates]

    def fillTemplate(self, tmpl):
        print('compiling template %s' % tmpl)
        with open(os.path.join(self.conf.model.templates, tmpl), 'r') as s:
            source = s.read()
            template = self.compiler.compile(source)
            outTmpl = template(self.data)
            tmplName = tmpl.split('.')[0]
            if tmplName == 'submitModel':
                ext = ".sh"
            else:
                ext = ".inp"
            outputFile = open(os.path.join(self.rundir, tmplName + ext), "w")
            [outputFile.write(str(line)) for line in outTmpl]
            outputFile.close()
