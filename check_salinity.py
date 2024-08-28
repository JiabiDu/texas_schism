#!/usr/bin/env python3
from pylib import *
close("all")
run='r262a3'
iday=60
itime=7
gd=loadz(f'../{run}/grid.npz').hgrid
#gd.plot()  #plot grid


#close('all')
#gd.plot()

#figure()
#gd.plot(fmt=1,value=gd.dp)


z=ReadNC(f'../{run}/outputs/temperature_{iday}.nc').temperature.val[itime,:,-1]
figure()
gd.plot(fmt=1,value=z,cmap='jet',clim=[15,25])
xlim([-94.76,-94.68])
ylim([29.755,29.795])
xlim([-95.4,-94.4])
ylim([29,30])
title(f'{run} day{iday:03d} itime{itime:02d}')
show()
