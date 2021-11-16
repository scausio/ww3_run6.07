import sys
import argparse
import os
import shutil
from glob import glob
from natsort import natsorted

def checkOutdir(path):
    if not os.path.exists(path):
        os.makedirs(path,0o755)
def main():
    outdir=os.path.join(path,'output')
    checkOutdir(outdir)
    dirs=natsorted([name for name in os.listdir(path) if os.path.isdir(path) if not name=='output'])
    for d in dirs:
        print ('working on %s'%d)
        run_path=os.path.join(path,d)
        fs=glob(os.path.join(run_path, 'ww3.*.nc'))
        if args.copy=='true':
            [shutil.copy(f,os.path.join(outdir,os.path.basename(f)))for f in fs]
        else:
            [shutil.move(f,os.path.join(outdir,os.path.basename(f)))for f in fs]

# this script move all netcdf produced by longrun in one folder.
parser = argparse.ArgumentParser()
parser.add_argument('--dir', '-d', help='directory')
parser.add_argument('--copy', '-c', help='copy files')

args = parser.parse_args()
print( args.dir)


#path='/work/opa/sc33616/ww3/blackSea_coup/3.1a/20180101'
path=args.dir
main() 
