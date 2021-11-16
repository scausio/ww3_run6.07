import os
import subprocess
import datetime
from configurationProvider import getConfigurationByID
import argparse


def longRun(startDate,duration,dayLimit,run ):
    runs=(int(duration)/int(dayLimit))
    startDate=datetime.datetime.strptime(startDate, "%Y%m%d")
    print(startDate)
    stopDate=startDate+ datetime.timedelta(int(duration))
    runningDate = startDate #+ datetime.timedelta(int(dayLimit)*int(run))
    deltaDays=(stopDate-runningDate).days
            
    if  deltaDays >=dayLimit:
            days=dayLimit
    else:
        days= deltaDays

    cmd = 'bsub -Is -P {p} -J longRun -o log.out -e log.err -q {queue} python main.py -s {date} -d {days} -r {run} -R {runs} -D {fullRun}'.format(queue=conf.model.serial_queue,date=runningDate.strftime( "%Y%m%d"),days=days,run=run,runs=runs,fullRun=duration,p=conf.model.project_queue)
    print(cmd)
    subprocess.call(cmd, shell=True)
            

parser = argparse.ArgumentParser()
parser.add_argument('--startDate', '-s', help='first run day')
parser.add_argument('--duration', '-d', help='run duration')
parser.add_argument('--run', '-r', help='run')
args = parser.parse_args()


startDate = args.startDate
duration = args.duration

if args.run:
    run = args.run
else:
    run = 0

conf = getConfigurationByID(os.path.dirname(os.getcwd()),'config')
longRun(startDate,duration,conf.model.output.dayLimit, run)
