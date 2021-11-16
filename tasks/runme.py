#!/usr/bin/python3
import argparse
from datetime import datetime, timedelta
import subprocess
from utils import getConfigurationByID, checkDir,setRunCMD
import os
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('--startDate', '-s', help='first run day')
parser.add_argument('--duration', '-d', help='run duration')
parser.add_argument('--endDate', '-e', help='last run day')
parser.add_argument('--dryRun', '-DR', help='dry run')
parser.add_argument('--workingDir', '-wd', help='run working dir')
args = parser.parse_args()

def get_runRange(args):
    startdate = args.startDate

    try:
        duration = int(args.duration)
        pass
    except:
        try:
            endDate = str(args.endDate)
            dt = datetime.strptime(endDate, "%Y%m%d") - datetime.strptime(startdate, "%Y%m%d")
            duration=dt.days
        except Exception as e:
            print (e, 'please add the end-date or duration of the run')
    return startdate, duration

def get_runRange_operational(args):
    startdate = args.startDate
    startdate = (datetime.strptime(startdate, "%Y%m%d")-timedelta(days=1)).strftime("%Y%m%d")
    try:
        duration = int(args.duration)+1
        pass
    except:
        try:
            endDate = str(args.endDate)
            dt = datetime.strptime(endDate, "%Y%m%d") - datetime.strptime(startdate, "%Y%m%d")
            duration=dt.days
        except Exception as e:
            print (e, 'please add the end-date or duration of the run')
    return startdate, duration



def setWorkingDir(config, startdate):
    conf_path = os.path.dirname(os.getcwd())
    conf = getConfigurationByID(conf_path, config)
    wd = os.path.join(conf.runpath, startdate)
    checkDir(wd)
    shutil.copy(os.path.join(conf_path, 'conf.yaml'), wd)
    checkDir(os.path.join(wd, 'output'))
    checkDir(os.path.join(wd, 'runs'))
    return conf,wd

def main():
    #NORMAL RUN
    startDate, duration = get_runRange(args)
    conf,wd = setWorkingDir('config',startDate)

    #OPERATIONAL RUN
    #startDate, duration = get_runRange_operational(args)
    #conf,wd = setWorkingDir('config',(datetime.strptime(startDate, "%Y%m%d")+timedelta(days=1)).strftime("%Y%m%d"))
    #
    cmd = setRunCMD(args)
    print (cmd.format(sla='',p=conf.model.project_queue,startdate=startDate, duration=duration, wd=wd, queue=conf.model.serial_queue))
    print ('dailyrun submitting')
    subprocess.call(cmd.format(sla='-Is',p=conf.model.project_queue,startdate=startDate, duration=duration, wd=wd, queue=conf.model.serial_queue), shell=True)

main()
