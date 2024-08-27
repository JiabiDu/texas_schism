#!/usr/bin/env python3
from pylib import *
close("all")

gd=loadz('grid.npz').hgrid
vd=loadz('grid.npz').vgrid
for itr, transect in enumerate(['../base/transect_channel.bp','../base/transect_shelf.bp','../base/transect_shelf_channel.bp','../base/transect_channel4seg.bp','../base/transect_channel4seg2.bp']):
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

    dis=[]
    for ist in arange(nst):
        dis.append(sqrt((bp.x[ist]-bp.x[max(ist-1,0)])**2+(bp.y[ist]-bp.y[max(ist-1,0)])**2)*111)
    dis=cumsum(dis)
    for ilay in arange(nlay):
        plot(dis,lzi[:,ilay]*-1,'-',color=[0,0,0,0.8],lw=0.5)

    for ist in arange(nst):
        plot(dis[ist]*ones(nlay),lzi[ist,:]*-1,'k-',lw=0.5)
    ylabel('Depth (m)')
    xlabel('Distance (km)')
    tight_layout()

    # add an inset to show the tramsect location
    if itr==0 or itr==3:
       axes(position=[0.2,0.25,0.3,0.4])
       gd.plot()
       gd.plot_bnd()
       xlim([min(bp.x)-0.03,max(bp.x)+0.03])
       ylim([min(bp.y)-0.03,max(bp.y)+0.03])
    else:
       axes(position=[0.5,0.25,0.4,0.5])
       gd.plot_bnd()
       #xlim([min(bp.x)-1,max(bp.x)+1])
       #ylim([min(bp.y)-1,max(bp.y)+1])
    plot(bp.x,bp.y,'r-')
    axis('off')
    savefig(f'{transect_name}_vgrid.png',dpi=300)
