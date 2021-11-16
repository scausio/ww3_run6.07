from datetime import datetime,timedelta
import subprocess
import numpy as np

start=20180101
days=365
step=100
spinup=3

runs=int(np.ceil(float(days)/step))

for i in range(runs):
    daysTorun=step
    counter=step*(i)
    print i
    if counter==0:
        startdate=(datetime.strptime(str(start),'%Y%m%d')+timedelta(counter)).strftime('%Y%m%d')
    else:
        startdate = (datetime.strptime(str(start), '%Y%m%d') + timedelta(counter-spinup)).strftime('%Y%m%d')
    if counter+step>days:
       
        daysTorun=days-counter + spinup
    cmd='bsub python main.py -s {stdate} -d {dur}'.format(stdate=startdate,dur=daysTorun)
    subprocess.call(cmd,shell=True)
    print cmd

