#!/usr/bin/env python3
import os
import sys
#-------------------------------------------------------------
#setup new folders based on previous run
#-------------------------------------------------------------
# sample call: ./copyrun.py RUN01a RUN01b

#input
base=sys.argv[1]
run=sys.argv[2]

#output dir
prj=os.getcwd().split('/')[-1] #project name, the last string cell in current working directory
opath=f'/scratch/group/galveston/jdu/{prj}/{run}/outputs'

def cpfile(sname,fname):
    if os.path.exists(fname): os.remove(fname) #need to remove thie file before linking
    os.symlink('../'+os.path.relpath(os.path.realpath(sname)),fname)

# remove run if existing
if os.path.exists(run+'/outputs/mirror.out'):
    print(run+' already exists and has model outputs')
    an=input(f'Are your sure to remove this {run} [y/n]?')
    if an not in ['Y','y']: sys.exit(0)

# copy run
os.system("rm -r {}; mkdir {}".format(run,run))
fnames=[i for i in os.listdir(base) if i.endswith(('in','nc','gr3','prop','ll','th','npz','sflux')) or i.startswith(('pschism',))]
for fname in fnames: #only use one level of sybmolic link
    cpfile(base+'/'+fname,run+'/'+fname)

# copy param.nml, run.*
os.system(f'cp {base}/param.nml {run}')
os.system(f'cp {base}/run* {run}')
os.system(f"cd {run}; rm -r outputs slurm-* *centers.bp screen.out")

# make output directory
os.system(f"cd {run}; mkdir -p  {opath}; ln -sf {opath} .")

# copy new input files; link depth = 1
if os.path.exists(f'Inputs/{run}'):
   fnames=[i for i in os.listdir(f'Inputs/{run}') if i.endswith(('in','nc','gr3','prop','ll','th','npz','sflux'))]
   for fname in fnames:
       cpfile(f'Inputs/{run}/'+fname,run+'/'+fname)

# copy flux.out, needed for hotstart
if os.path.exists(f'../{base}/outputs/flux.out'):
    os.system(f"cd {run}; cp ../{base}/outputs/flux.out outputs/") #needed for hotstart

print(f'Succesfully copy {base} to {run} for project {prj}')
sys.exit()
