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
host='grace'
#if len(sys.argv)>3: host=sys.argv[3]
print(base,run)

#output dir
#opath='/sciclone/scr10/jiabi/GoM'
#if host in ['frontera','f','front']: opath='/scratch1/08709/jiabi/GoM'  #to run on frontera

#copy run
os.system("rm -rf {}; mkdir {}".format(run,run));
#os.system("cd {}; cp -a ../{}/* ./".format(run,base))
os.system(f"cd {run}; rsync -av --exclude={'outputs'} ../{base}/* .")
#odir="{}/{}".format(opath,run);
os.system(f"cd {run}; rm -r outputs slurm-* *centers.bp screen.out; mkdir outputs")
print(f'Succesfully copy {base} to {run} on {host}')
sys.exit()


