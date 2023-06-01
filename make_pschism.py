#!/usr/bin/env python3
import os
import sys
cmd='cd ~/schism/mk; ln -sf Make.defs.grace Make.defs.local; cd ../src; make clean; make pschism'
print(cmd); os.system(cmd)

#add git tag
tag=os.popen('cd ~/schism; git log').read().split('\n')[0].split()[1][:8]

#get the host name
host=os.popen('echo $HOST').read().split()[0].upper()

#add git tag and host name to the ptrack3.WW
fname=[i for i in os.listdir('/home/jdu/schism/src/') if i.startswith('pschism_GRACE_Intel')]
if len(fname)>0: 
    exe=fname[0]
else:
    sys.exit('not compile successfully')

cmd=f'mv ~/schism/src/{exe} ./{exe}.{tag}'
print(cmd); os.system(cmd)
