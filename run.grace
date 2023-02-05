#!/usr/bin/env python3
'''
  atuo script to submit SCHISM batch jobs on sciclone/james
  note: environmental variables "run_schism" is used to load correct modules 
'''
import os
import sys
from pylib import get_hpc_command
#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
#code='pschism_FRONTERA_ICM_OLDIO_PREC_EVAP_TVD-VL.3ff3e18d'
#code='pschism_VORT_ICM_OLDIO_PREC_EVAP_TVD-VL.3754530e'
#code='pschism_FRONTERA_OLDIO_TVD-VL.07e4b6e3'
code='pschism_GRACE_Intel_VL.cad100e5 4'
#code='pschism_BORA_OLDIO_TVD-VL.07e4b6e3'
runs=[]  #used to submit multiple jobs

#resource requst 
walltime='16:00:00'
qnode='grace'; nnode=6; ppn=48; mem='32G' 
#qnode='x5672'; nnode=2; ppn=8      #hurricane, ppn=8
#qnode='bora'; nnode=1; ppn=20      #bora, ppn=20
#qnode='vortex'; nnode=20; ppn=12    #vortex, ppn=12
#qnode='femto'; nnode=6; ppn=32     #femto,ppn=32
#qnode='potomac'; nnode=4; ppn=8    #ches, ppn=12
#qnode='james'; nnode=5; ppn=20     #james, ppn=20
#qnode='skylake'; nnode=2; ppn=36   #viz3,skylake, ppn=36
#qnode='haswell'; nnode=2; ppn=2    #viz3,haswell, ppn=24,or 28
#qnode='frontera'; nnode=10; ppn=56  #frontera, ppn=56 
#qnode='mistral'; nnode=5; ppn=36    #mistral, ppn=36
qname='flex'                        #partition name (needed for frontera)

scrout='screen.out'; bdir=os.path.abspath(os.path.curdir); jname=os.path.basename(bdir)
#-----------------------------------------------------------------------------
#submit jobs first (qsub)
#-----------------------------------------------------------------------------
if os.getenv('run_schism')==None: 
   scode=get_hpc_command(sys.argv[0],bdir,jname,qnode,nnode,ppn,walltime,scrout,0,'run_schism',qname=qname,mem=mem)
   print(scode); os.system(scode); os._exit(0)

#-----------------------------------------------------------------------------
#running parallel jobs (mpirun) 
#-----------------------------------------------------------------------------
if os.getenv('run_schism')!=None: 
   fmt=1; bdir=os.getenv('run_schism').split()[0]; os.chdir(bdir)
   scode=get_hpc_command(code,bdir,jname,qnode,nnode,ppn,walltime,scrout,1,'run_schism',mem=mem)

   if len(runs)==0: #single jobs 
      print(scode); sys.stdout.flush(); os.system(scode)
   else: #multiple jobs: put script outside of run dir
      for m,run in enumerate(runs): 
          print(scode); os.chdir(bdir+'/'+run); os.system(scode)
   sys.exit(0) if qnode in ['bora'] else os._exit(0)
