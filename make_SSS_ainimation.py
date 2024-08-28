#!/usr/bin/env python3
from pylib import *
from pcom import *
close("all")
run='run19m3'
reftime=datenum(2018,1,1)
varname='salinity'; clims=[12,37]
replot=True
#varname='temperature'; clims=[20,30]

os.makedirs(f'{run}/figs',exist_ok=True)
gd=read_schism_hgrid(f'../{run}/hgrid.ll')
# figure(figsize=[11,4])
# for ifile in [6]:
#     C=ReadNC(f'../{run}/outputs/{varname}_{ifile}.nc',fmt=1) 
#     mtime=C['time'][:]
#     for itime in arange(0,len(mtime)):
#         figname='{}/figs/SSS_{:06.2f}.png'.format(run,mtime[itime]/86400)
#         if fexist(figname): print(figname,'exists'); continue
#         clf()
#         axes(position=[0.05,0.12,0.55,0.85])
#         gd.plot(fmt=1,value=C[varname][itime,:,-1],cmap='gist_rainbow_r',clim=clims,cb=None);
#         xlabel('Lon (deg)')
#         ylabel("Lat (deg)")
#         rtext(0.05,0.95,'{:.2f} days since 2018-01-01'.format(mtime[itime]/86400))
        
#         axes(position=gap(rel=[1.1,0,0.6,1]))
#         gd.plot(fmt=1,value=C[varname][itime,:,-1],cmap='gist_rainbow_r',clim=clims,ticks=6);
#         xlim([-95.2,-94.4]); ylim([29,29.8])
#         print(figname);savefig(figname,dpi=200)
#         #show(); pp
# make_gif(image_dir=f'{run}/figs/',begin_string='SSS_',gif=f'{run}/SSS.gif',duration=100)

add_timeseries=True
if add_timeseries:
    reftime=datenum(2018,1,1)
    M=loadz(f'{run}/salt.npz')
    fp=M.bp.station=='MIDG'
    fp0=nonzero(fp)[0]
    tmpz=M.bp.z[fp0]
    tmp=abs(tmpz-3) #find the depth closest to the needed depth. 
    tfp=fp0[tmp==tmp.min()]
    olon,olat=M.bp.x[tfp],M.bp.y[tfp]

    O=loadz('../Observations/twdb_wq/npz/twdb_salinity.npz')
    mx=M.time+reftime
    my=M.salinity[tfp].ravel()
    my=lpfilt(M.salinity[tfp].ravel(),median(diff(M.time)),0.48) #0.48==50hour
    
    fp=(O.station=='MIDG')*(O.var=='seawater_salinity')*(O.time>M.time.min()+reftime)*(O.time<M.time.max()+reftime)
    ox,oy=O.time[fp],O.data[fp]
    oy=lpfilt(O.data[fp],median(diff(O.time[fp])),0.48)

figure(figsize=[11,4])
for ifile in arange(1,37):
    C=ReadNC(f'../{run}/outputs/{varname}_{ifile}.nc',fmt=1)
    mtime=C['time'][:]
    for itime in arange(11,len(mtime),12):
        figname='{}/figs/SSS2_{:06.2f}.png'.format(run,mtime[itime]/86400)
        if fexist(figname) and not replot: print(figname,'exists'); continue
        clf()
        axes(position=[0.05,0.12,0.55,0.85])
        gd.plot(fmt=1,value=C[varname][itime,:,-1],cmap='gist_rainbow_r',clim=clims,cb=None);
        # add the shelf break
        xlabel('Lon (deg)')
        ylabel("Lat (deg)")
        rtext(0.05,0.95,'{:06.2f} days since 2018-01-01'.format(mtime[itime]/86400))
         

        axes(position=gap(rel=[1.1,0,0.6,1]))
        gd.plot(fmt=1,value=C[varname][itime,:,-1],cmap='gist_rainbow_r',clim=clims,ticks=6);
        xlim([-95.2,-94.4]); ylim([29,29.8])
        plot(olon,olat,'k.',ms=10)
        text(olon-0.05,olat+0.02,'MIDG')


        if add_timeseries:
            axes(position=[0.3,0.2,0.25,0.2])
            plot(ox-reftime,oy,'x',color=[0,0,0,0.5],lw=0.5,ms=1)
            fp=mx-reftime<=mtime[itime]/86400
            plot(mx[fp]-reftime,my[fp],'-',color=[1,1,1,1],lw=1)
            xlim([0,365])
            ylim([0,27])
            gca().patch.set_alpha(0.0)
            text(150,10,'Obs',color='k',fontweight='bold')
            text(150,5,'Model',color='w',fontweight='bold')
            title('Salinity at MIDG',color='w')
            gca().spines['top'].set_visible(False)
            gca().spines['right'].set_visible(False)
        print(figname);savefig(figname,dpi=150)
        #show(); pp

make_gif(image_dir=f'{run}/figs/',begin_string='SSS2_',gif=f'{run}/SSS2.gif',duration=100)
