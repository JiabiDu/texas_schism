#!/usr/bin/env python3
from pylib import *
close('all')
#!! remember to add rivers.bp after using this script, which will be used for source_sink
#----------------------------------------------------------------------
#inputs
#----------------------------------------------------------------------
#grd=[i for i in os.listdir('.') if i.endswith('.2dm')][0]
grd=os.getcwd().split('/')[-1]+'.2dm'
min_depth=1
#dp_upper_atcha=10
river_bp='../base/rivers_V2.bp'
if not fexist(grd): sys.exit(f'{grd} not exist; file name must be consistent with the folder name')
prj='epsg:4326'  #projection of grd
bxy='../base/bxy3.bp' #lon&lat of bnd nodes
headers=("etopo1","crm_3arcs","cdem13_",'ncei19_2019v1_',"cb_bay_dem_v3.1_ll") #DEM headers
headers=("etopo1","ncei19_MS_LA", "ncei19_AL_nwFL","ncei19_TX")
bdir=r'/scratch/user/jdu/DEM/npz'  #directory of DEM data
ogd=read_schism_hgrid('../v000/hgrid.gr3') #bathymetry in east Miss is used.  
use_lsc2=False

#----------------------------------------------------------------------
#get hgrid.gr3, hgrid.ll and vgrid.in
#----------------------------------------------------------------------
#read 2dm, do projection, compute bndfile, and split bad quads
print(f'converting {grd} to *.gr3')
gd=sms2grd(grd); lx,ly=gd.proj(prj,fmt=1); px,py=gd.proj(prj,'epsg:26915',fmt=1)
gd.check_skew_elems(angle_min=15,fname='skew_elem_15.bp')
gd.check_skew_elems(angle_min=10,fname='skew_elem_10.bp')
gd.check_skew_elems(angle_min=5,fname='skew_elem_5.bp')
gd.check_skew_elems(angle_min=2,fname='skew_elem_2.bp')
gd.x,gd.y=px,py; gd.split_quads(angle_min=60,angle_max=120) #split bad quads
gd.x,gd.y=lx,ly; gd.compute_bnd(bxy) #define boundary


#load bathymetry
print('loading bathymetry to grid')
for header in headers: 
    DEMs=array([i for i in os.listdir(bdir) if i.startswith(header)])
    for n,DEM in enumerate(DEMs): 
        dp,ip=load_bathymetry(gd.x,gd.y,'{}/{}'.format(bdir,DEM),fmt=1); gd.dp[ip]=-dp

#minimum depth in regions, etc.  (todo: regions)
print(f'set min depth to {min_depth}m')
gd.dp[gd.dp<min_depth]=min_depth

# change grid bathymetry east to -93.6 based on original one
print('use old grid bathymetry for region with x>-93.6')
fp=near_pts(c_[gd.x,gd.y],c_[ogd.x,ogd.y])
newdp=ogd.dp[fp]
fp=gd.x>-93.6
gd.dp[fp]=newdp[fp]

#set a larger depth for the upper Atchafalaya River
#bp=read_schism_reg('upper_atchafalaya.reg')
#fp=nonzero(inside_polygon(c_[gd.x,gd.y],bp.x,bp.y)==1)[0]
#for ip in fp:
#    gd.dp[ip]=max(dp_upper_atcha,gd.dp[ip])

#gd.dp[fp]=concatenate((gd.dp[fp],ones(sum(fp))*dp_upper_atcha))

#save grid
gd.save('hgrid.gr3',fmt=1)

if not use_lsc2:
    os.system('cp ../base/vgrid.in .')
else:
    print('generating vgrid.in')
    code='ifort -Bstatic -O2 -CB -o gen_vqs gen_vqs.f90 schism_geometry.f90'
    os.system('ln -sf ../base/*.f90 .; {}; ./gen_vqs; cp vgrid.in vgrid.in.old'.format(code)) #gen vgrid.in
    read_schism_vgrid('vgrid.in').write_vgrid() #change vgrid format

#clean folders
os.system('ln -sf hgrid.gr3 hgrid.ll')
os.system('rm estuary.gr3 fort.* gen_vqs* schism_geometry* vgrid_master_*')
save_schism_grid()

os.system(f'ln -sf {river_bp} rivers.bp')

# check the grid
gd=loadz('grid.npz').hgrid
figure(figsize=[10,5])
gd.plot(fmt=1,clim=[0,20]); gd.plot(); plot(gd.x[gd.iobn[0]],gd.y[gd.iobn[0]],'bo');
xlim([gd.x.min(),gd.x.max()])
ylim([gd.y.min(),gd.y.max()])
bp=read_schism_bpfile('rivers.bp')
plot(bp.x,bp.y,'ro')
tight_layout()
savefig('hgrid.png',dpi=300)
show()
