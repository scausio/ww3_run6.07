# install conda environment using environment.yml in ww3_run6.07/setup

conda env create -f environment.yml

# load libraries and environment

- set properly the PATH variable in setup.sh to point to your ww3 directories

-source ww3_run6.07/setup/setup.sh

# set the run

- add your .msh file in ww3_run6.07/data/mesh
- for nested run, add your list of boundary nodes file in ww3_run6.07/data/boundary
- fill properly the conf.yaml file

# run

cd to tasks folder

bsub -o log.out -e log.err  python runme.py -s [startdate] -d [days] or -e [enddate]
