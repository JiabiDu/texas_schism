#!/usr/bin/env python3
from pylib import *
close("all")
run='run19a'
ifile=10
xlims=[-94.73-0.2,-94.73+0.3]
ylims=[29.32-0.2,29.32+0.2]
xlims=[-94.74,-94.64]
ylims=[29.29,29.37]
gd=loadz(f'../{run}/grid.npz').hgrid
U=ReadNC(f'../{run}/outputs/horizontalVelX_{ifile}.nc',fmt=1)
V=ReadNC(f'../{run}/outputs/horizontalVelY_{ifile}.nc',fmt=1)
itime=30
fp=(gd.x>xlims[0])*(gd.x<xlims[1])*(gd.y>ylims[0])*(gd.y<ylims[1])
u=U['horizontalVelX'][itime,fp,-1]
v=V['horizontalVelY'][itime,fp,-1]
speed=sqrt(u**2+v**2)
quiver(gd.x[fp], gd.y[fp], u, v, speed,width=0.002,pivot='tip',cmap='jet',scale=20)
gd.plot_bnd()
xlim(xlims); ylim(ylims)
clim([0,1.5])
colorbar()
savefig('tmp/flowGalveston.png',dpi=300)
show()

