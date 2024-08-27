#!/usr/bin/env python3
from pylib import *
from glob import glob
close('all')
#!! remember to add rivers.bp after using this script, which will be used for source_sink
#----------------------------------------------------------------------
#inputs
#----------------------------------------------------------------------
#grd=[i for i in os.listdir('.') if i.endswith('.2dm')][0]
if len(glob('*.2dm'))>1: sys.exit('there are more than one 2dm file')
grd=glob('*.2dm')[0]
min_depth=1
river_bp='../base/rivers_V2.bp'
if not fexist(grd): sys.exit(f'{grd} not exist; file name must be consistent with the folder name')
prj='epsg:4326'  #projection of grd
bxy='../base/bxy4.bp' #lon&lat of bnd nodes
headers=("etopo1","crm_3arcs","cdem13_",'ncei19_2019v1_',"cb_bay_dem_v3.1_ll") #DEM headers
headers=("etopo1","ncei19_LA_MS", "ncei19_AL_nwFL","ncei19_TX",'TWDB2','CESWG')
bdir=r'/scratch/user/jdu/DEM/npz'  #directory of DEM data
use_lsc2=True
vgrid='../vgrid_v5' #vgrid version
one_layer_trinity=True

#----------------------------------------------------------------------
#get hgrid.gr3, hgrid.ll and vgrid.in
#----------------------------------------------------------------------
#read 2dm, do projection, compute bndfile, and split bad quads
if not fexist('hgrid.gr3'):
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
        DEMs.sort()
        print('loading DEM',header)
        for n,DEM in enumerate(DEMs): 
            dp,ip=load_bathymetry(gd.x,gd.y,'{}/{}'.format(bdir,DEM),fmt=1); gd.dp[ip]=-dp
    
    #minimum depth in regions, etc.  (todo: regions)
    print(f'set min depth to {min_depth}m')
    gd.dp[gd.dp<min_depth]=min_depth
    
    #set min depth of 2m for trinity channel (to resolve abnormal speed)
    #bp=read_schism_reg('../base/trinity_channels.reg')
    #fp=(inside_polygon(c_[gd.x,gd.y],bp.x,bp.y)==1)*(gd.dp<=2)
    #gd.dp[fp]=2
    
    bp=read_schism_reg('../base/upperGalveston01.reg') #some abnoraml high after including TWDB2
    fp=(inside_polygon(c_[gd.x,gd.y],bp.x,bp.y)==1)*(gd.dp>5)
    gd.dp[fp]=1.2
    print(f'number of abnormal deep point to be corrected: {sum(fp)}')
        
    #set minimum dpeth at the open boundary to 5
    for ii in gd.iobn[0]:
        gd.dp[ii]=max(5,gd.dp[ii])
    
    #save grid
    gd.save('hgrid.gr3',fmt=1)
else:
    print('hgrid.gr3 already exists')

# vertical grid
if not use_lsc2:
    os.system('cp ../base/vgrid.in .')
else:
    gd=read_schism_hgrid('hgrid.gr3'); dp0=gd.dp.copy()
    bp=read_schism_reg(vgrid+'/estuary.reg')
    gd.dp[:]=0; gd.dp[inside_polygon(c_[gd.x,gd.y],bp.x,bp.y)==1]=1; gd.save('estuary.gr3')
    gd.dp=dp0 #put back the original depth 
  
    print('generating vgrid.in')
    cmd=f'ln -sf {vgrid}/*.f90 .; ifort -Bstatic -O2 -CB -o gen_vqs gen_vqs.f90 schism_geometry.f90; ./gen_vqs; cp vgrid.in vgrid.in.old'
    print(cmd); os.system(cmd)
  
    vd=read_schism_vgrid('vgrid.in')
    print('make vertical grid evenly distributed when dp<5')
    fp=nonzero(gd.dp<5)[0]
    for ind in fp:
        sigma=vd.sigma[ind].copy()
        nlevel=sum(sigma>-1.0)+1 #number of levels
        sigma[-nlevel:]=linspace(-1,0,nlevel)
        vd.sigma[ind]=sigma 
    
    # use one layer in the trinity
    if one_layer_trinity:
        bp=read_schism_reg('../base/trinity_mouth.reg')
        fp=nonzero(inside_polygon(c_[gd.x,gd.y],bp.x,bp.y))[0] #grids inside the reg
        print(f'set two levels (one layer) for {len(fp)} nodes in Trinity')
        for ind in fp:
            vd.sigma[ind][:]=-1.0
            vd.sigma[ind][-2:]=linspace(-1.0,0,2)

    vd.write_vgrid() #change vgrid format
    os.system('rm estuary.gr3 fort.* gen_vqs *.f90 *.mod vgrid_master*') #clean folder

#make a symbolic link and save grid into grid.npz
os.system('ln -sf hgrid.gr3 hgrid.ll')
save_schism_grid()
gd=loadz('grid.npz').hgrid
vd=loadz('grid.npz').vgrid

# get rivers.bp
fnames=[i for i in os.listdir('../rivers/') if i.endswith('.bp') and i!='rivers.bp']
x,y,river_name=[],[],[]
mindist=[]
for fname in fnames:
    bp=read_schism_bpfile('../rivers/'+fname)
    x.append(bp.x[0])
    y.append(bp.y[0])
    ig=near_pts(c_[bp.x[0],bp.y[0]],c_[gd.x,gd.y])
    mindist.append(sqrt((gd.x[ig]-x[-1])**2*cos(y[-1]*pi/180)**2+(gd.y[ig]-y[-1])**2)*11e3)
    river_name.append(' '.join(fname.split('_')[0:-1]))
x,y,river_name,mindist=array(x),array(y),array(river_name),array(mindist).ravel()
# find those within the model domain or has a min distance of less than 200m from domain
ie,ip,acor=gd.compute_acor(c_[x,y])
fp=(ie!=-1)|(mindist<200)
x,y,river_name=x[fp],y[fp],river_name[fp]

#sort by longitude
fp=argsort(x)
x,y,river_name=x[fp],y[fp],river_name[fp]
print('writing rivers.bp')
bp=schism_bpfile()
bp.nsta=len(x);bp.x=x;bp.y=y;bp.z=zeros(len(x));bp.station=river_name
bp.write_bpfile('rivers.bp')

# plot the grid and show the river locations 
figure(figsize=[10,5])
gd.plot(fmt=1,clim=[0,20]);
plot(gd.x[gd.iobn[0]],gd.y[gd.iobn[0]],'bo'); #open boundary
xlim([gd.x.min(),gd.x.max()])
ylim([gd.y.min(),gd.y.max()])
plot(bp.x,bp.y,'ro',ms=2) #rivers
for m,(xi,yi) in enumerate(zip(bp.x,bp.y)):
    text(xi,yi,f'{m+1}',color='r')
tight_layout()
savefig('hgrid.png',dpi=300)

#%% plot the vertical grid across multiple transects
for itr, transect in enumerate(['../base/transect_shelf.bp','../base/transect_shelf_channel.bp','../base/transect_channel.bp','../base/transect_channel4seg.bp','../base/transect_channel4seg2.bp','../base/transect_trinity.bp']):
    transect_name=transect.split('/')[-1].split('.')[0]
    bp=read_schism_bpfile(transect)
    bind=near_pts(c_[bp.x,bp.y],c_[gd.x,gd.y])
    if vd.ivcor==2:
         lzi=abs(compute_zcor(vd.sigma,gd.dp[bind],ivcor=2,vd=vd))
    else:
         lzi=abs(compute_zcor(vd.sigma[bind],gd.dp[bind]));

    figure(figsize=(6,3))
    nst=len(bp.x)
    nlay=lzi.shape[1]

    dis=cumsum([sqrt((bp.x[ist]-bp.x[max(ist-1,0)])**2+(bp.y[ist]-bp.y[max(ist-1,0)])**2)*111 for ist in range(nst)])
    for ilay in arange(nlay):
        plot(dis,lzi[:,ilay]*-1,'-',color=[0,0,0,0.8],lw=0.5)

    for ist in arange(nst):
        plot(dis[ist]*ones(nlay),lzi[ist,:]*-1,'k-',lw=0.5)
    ylabel('Depth (m)'); xlabel('Distance (km)')
    tight_layout()

    # add an inset to show the tramsect location
    if itr>=2:
       axes(position=[0.15,0.25,0.3,0.4])
       gd.plot()
       gd.plot_bnd()
       xlim([min(bp.x)-0.01,max(bp.x)+0.01])
       ylim([min(bp.y)-0.0075,max(bp.y)+0.0075])
    else:
       axes(position=[0.5,0.25,0.4,0.5])
       gd.plot_bnd()
    plot(bp.x,bp.y,'k.',lw=0.5,ms=5)
    axis('off')
    savefig(f'{transect_name}_vgrid.png',dpi=300)
print('completed!')
show()
