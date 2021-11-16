import os
import subprocess
from glob import glob
from preprocessing import Preprocessing
from utils import Paths, Checks, setRunCMD
from datetime import datetime,timedelta
from utils import getConfigurationByID, checkDir, nextDay,checkFunctionCompleted,checkDayCompleted,dataFormatter
import argparse
from postprocessing import Postprocessing,Regridding
import shutil
from bjobs import Bjobs

parser = argparse.ArgumentParser()
parser.add_argument('--startDate', '-s', help='first run day')
parser.add_argument('--duration', '-d', help='run duration')
parser.add_argument('--dryRun', '-DR', help='dry run')
parser.add_argument('--workingDir', '-wd', help='run working dir')
args = parser.parse_args()


class NextRun():
    def __init__(self,conf,totalDuration,runDuration, paths,args,submitter):
        self.conf=conf
        self.submitter=submitter
        self.totalDuration=totalDuration
        self.runDuration=runDuration
        self.workingdir=paths.workingDir
        self.paths=paths

        if self.totalDuration>self.runDuration:
            self.buildCMD(args)
            self.lastRun=False
        else:
            self.lastRun=True
            self.cmd='run over'

    def buildCMD(self, args):
        nextDuration = self.totalDuration - self.runDuration
        nextStartDay = dataFormatter(args.startDate)+timedelta(days=self.runDuration)

        wd = args.workingDir
        cmd= setRunCMD(args)
        self.cmd = cmd.format(sla='' ,p=self.conf.model.project_queue,queue=self.conf.model.serial_queue, startdate=nextStartDay.strftime('%Y%m%d'), duration=nextDuration, wd=wd)
        return self.cmd
    def run(self):
        if self.totalDuration >self.runDuration:
            os.chdir(os.path.join(self.conf.base, 'tasks'))
            self.submitter.systemCommand(self.cmd)
        else:
            print (self.totalDuration, self.runDuration, self.workingdir)
            #os.remove(glob(os.path.join(self.workingdir,'restart*.ww3'))[-1])
            if self.conf.post.regrid.flag.upper() in ['T','TRUE','Y','YES']:
                regridder=Regridding(self.conf,self.paths)
                for area in self.conf.post.regrid.area:
                    print(area)
                    vars = self.conf.post.regrid.area[area].variables
                    if vars:
                        regridder.regriddingMain(area)
            print ('run over')

class ModelRun():
    def __init__(self, paths, checks,submitter,nextBsub):
        self.paths=paths
        self.submit(checks,submitter,nextBsub)

    def submit(self,  checks, submitter,nextBsub):
        if str(args.dryRun).upper() in ['Y','YES']:
            print('dry run activated')
        else:
            checkFunctionCompleted(checks.SHELcomplete, submitter.model)
        f=open(os.path.join(self.paths.runDir, 'bsub.txt'), 'w')
        f.write("{}".format(nextBsub))
        f.close()

def  setRunDirectories(paths):
    print ('setting directories')
    [checkDir(getattr(paths, directory)) for directory in paths.__dict__.keys()]
    os.chdir(paths.runDir)
    return paths.runDir

class Submit():
    def __init__(self,conf,paths):
        self.paths=paths
        self.conf=conf
    def model(self):
        os.chdir(self.paths.runDir)
        self.systemCommand('chmod +x submitModel.sh')
        print ('Submitting model ...')
        Bjobs(['bsub < submitModel.sh'] ,self.conf.model.project_queue)#,check=True,shell=True)
        #cmd='bsub -K -q {queue} -n {n} mpirun {exe}/ww3_shel'.format(queue=self.conf.model.multiprocess.queue,n=self.conf.model.multiprocess.cores,exe=self.conf.model.executable)
        #subprocess.run(cmd,check=True,shell=True)
    def systemCommand(self,cmd):
        print ('Submitting %s'%cmd)
        subprocess.run(str(cmd), check=True,shell=True)


def main():
    startDate=str(args.startDate)
    duration=int(args.duration)
    conf = getConfigurationByID(args.workingDir, 'config')
    if duration> conf.model.limitDuration:
        runDuration=conf.model.limitDuration
    else:
        runDuration=duration

    paths= Paths(conf, startDate, args)
    submitter = Submit(conf,paths)
    nextrun = NextRun(conf,duration,runDuration,paths, args, submitter)
    if not checkDayCompleted(paths,startDate):
        runDir = setRunDirectories(paths)
        checks = Checks(runDir)
        os.chdir(paths.runDir)
        Preprocessing(conf, runDir,runDuration, args, checks, submitter,nextrun.lastRun).main()
        ModelRun(paths,checks,submitter,nextrun.cmd)
        Postprocessing(conf,paths,args,checks,submitter)
    nextrun.run()

main()
