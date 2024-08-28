#!/usr/bin/env python3
from pylib import *
from pcom import *
close("all")

grid='/scratch/user/jdu/GoM/Grids/V2.7.0/'
gd=loadz(f'{grid}/grid.npz').hgrid
bp=read_schism_bpfile(f'{grid}/rivers.bp')

figure(figsize=[8,5])
gd.plot()
plot(bp.x,bp.y,'ro',ms=2) #rivers
xlabel('Lon (deg)')
ylabel('Lat (deg)')
xlim([gd.x.min(),gd.x.max()])
ylim([gd.y.min(),gd.y.max()+0.1])
tight_layout()
adj_figsize()
show()
savefig('final_report/grid.png',dpi=600)
