#!/usr/bin/env python3
from pylib import *   
from pcom import *    
from fgom import *
close("all") #close all existing figure windows
os.system('clear')  #clean the screen
run='r270b'
brun=None
reftime=datenum(2018,1,1)
show_location=True;

p=zdata()  #a container class, used to install key parameters
p.salt_plot         = 1 
p.elev_plot         = 1
p.current_plot      = 1
p.shelfcurrent_plot = 1
p.temp_plot         = 1
p.harmonic          = 0
p.tidal_names=np.array(['K1','O1','Q1','P1','M2','S2','N2','K2'])
p.salt_obs='../Observations/twdb_wq/npz/twdb_salinity.npz'
p.salt_nerrs='../Observations/nerr/NERRS.npz'
p.elev_obs='../Observations/noaa_wl_navd/npz/hourly_height_2000_2022.npz'
p.elev_obs_info='../Observations/noaa_wl_navd/npz/stainfo.npz'
p.temp_obs='../Observations/noaa_wl_temp/npz/water_temperature_2000_2022.npz'

gd=read_schism_hgrid(f'../{run}/hgrid.ll')
os.makedirs(run+'/figs',exist_ok=True)
ioff()

#%%
if p.salt_plot:
    M=loadz(run+'/salt.npz')
    B=None
    if brun!=None and type(brun)!=list: B=loadz(brun+'/salt.npz')
    if type(brun)==list: B=[loadz(brun+'/salt.npz') for brun in brun] #if multiple base runs to compare
    
    # plot salinity for different bays
    O=loadz(p.salt_obs)
    O.var[O.var=='seawater_salinity']='salinity'
    bay_station={'Galveston':[('BOLI','MIDG','FISH','TRIN'),(3,3,1.5,1.5),None],
                 'Sabine':[('SAB1','SAB2'),(2,2),None],
                 'Laguna_Corpus':[['LMA2','BIRD','SPCG','INPT'],(1,1,1,1),[10,40]],
                 'Saint':[('GEA-S','GEA-B','DELT','MOSQ'),(1,1,1,1,1,1),None],
                 'Matagorda':[('TPAL','6985','6990','6980','EMB'),(1,1,1,1,1),None],} 
    for fmt in [1,0]: #full and subtidal
        for bay in bay_station:
            plot_salinity_one_layer(M,B,O,fmt=fmt,stations=bay_station[bay][0],deps=bay_station[bay][1],bay=bay, run=run,brun=brun,reftime=reftime, gd=gd, 
                                    show_location=show_location,mlg=-1,ylims=bay_station[bay][2])
            pppp

    # salinity at NERRS stations   
    stations=['CE','CW','MB','AB']; deps=(1.5,1.5,1.5,1.5)
    O=loadz(p.salt_nerrs)
    O.var=O.var.astype('<U8'); O.var[O.var=='Sal']='salinity'
    O.lon={i:O.lon[O.station==i].mean()*-1 for i in stations} # a dictionary storing the lon,lat for the stations
    O.lat={i:O.lat[O.station==i].mean() for i in stations}
    plot_salinity_one_layer(M,B,O,fmt=1,stations=stations,deps=deps,bay='Mission_Aransas',run=run,brun=brun,reftime=reftime, gd=gd,show_location=show_location,mlg=0,ylims=None)

#%% fmt=0,1,2 for full, subtidal, tidel
if p.elev_plot:
    M=loadz(run+'/elevation.npz')
    B=None
    if brun!=None: B=loadz(brun+'/elevation.npz')
    O=loadz(p.elev_obs)

    # get station information
    E=loadz(p.elev_obs_info)
    lon={i:j for i,j in zip(E.station,E.lon)}
    lat={i:j for i,j in zip(E.station,E.lat)}
    station_name={i:j for i,j in zip(E.station,E.station_name)}

    #shoreline stations, full, tidal, and subtidal
    stations=['8779749','8775870','8775241','8773767','8772471','8770613','8771341','8770822'];
    regname='shore'
    for signal in ['full','tidal','subtidal']:
        if signal in ['full','tidal']: xlims=xlims=[datenum(2018,2,1),datenum(2018,2,20)]
        if signal=='subtidal':xlims=None
        plot_elevation(stations,M,B,O,row=4,col=2,margin=[0.1,0.02,0.06,0.1],dxy=[0.08,0.00],signal=signal,
                       xlims=xlims,ylims=[-0.5,0.5],gd=gd,
                    run=run,brun=brun,reftime=reftime,station_name=station_name,lon=lon,lat=lat,
                    regname=regname,paper=True,adj=False) #full signal

    #all stations, full, subtidal and tidal 
    stations=[i for i in M.bp.station if i.startswith('87')]
    regname=''
    for singal in ['full','tidal','subtidal']:
        if signal in ['full','tidal']: xlims=xlims=[datenum(2018,2,1),datenum(2018,2,20)]
        if signal=='subtidal':xlims=None
        stations=['8779749','8775870','8775241','8773767','8772471','8770613','8771341','8770822'];
        plot_elevation(stations,M,B,O,row=4,col=2,margin=[0.1,0.02,0.06,0.1],dxy=[0.08,0.00],signal=signal,
                       xlims=xlims,ylims=[-0.5,0.5],gd=gd,
                    run=run,brun=brun,reftime=reftime,station_name=station_name,lon=lon,lat=lat,
                    regname='shore',paper=True,adj=False)
    
#%%
if p.current_plot:
    O=loadz('../Observations/noaa_current/noaa_current_gom_2018.npz')
    M=loadz(run+'/current.npz')
    M.speed=sqrt(M.horizontalVelX**2+M.horizontalVelY**2)
    B=None
    if brun!=None: 
        B=loadz(brun+'/current.npz')
        B.speed=sqrt(B.horizontalVelX**2+B.horizontalVelY**2)
    gd=read_schism_hgrid(f'../{run}/hgrid.ll')
    stations=['STX1819','STX1801','STX1802','g06010','mg0201','sn0301','sn0501','sn0701']
    for fmt in [0,1,2]: 
        plot_current(M,O,B,p=p,zoom=5,stations=stations,col=1,row=4,reg='shore',fmt=fmt,
                     run=run,brun=brun,gd=gd,reftime=reftime) #zoom in view the first obs-available 30 days

    #plot_current(M,O,B,p=p,zoom=30) #zoom in view the first obs-available 30 days
    #plot_obs_station(M,gd=gd,var='current',p=p)
    #plot_current_direction_hist(M,O,p=p)
    #plot_current_speed_hist(M,O,p=p)
    #plot_current(M,O,B,p=p)
#%%
if p.shelfcurrent_plot:
    O=loadz('../Observations/ndbc_current/ndbc_current.npz')
    M=loadz(run+'/current.npz')
    gd=loadz(f'../{run}/grid.npz').hgrid
    B=None
    if brun!=None: B=loadz(brun+'/current.npz')
    for subtidal in [False,True]:
        xlims=[M.time.min(),M.time.max()]
        #if subtidal==False: xlims=[100,200]
        plot_shelf_current(O,M,B,run=run,brun=brun,gd=gd,stations=[42043,42050],station_nickname=['B','F'],xlims=xlims,
                       fmt=0,show_location=True,subtidal=subtidal)

#%% 
if p.temp_plot:
    O=loadz(p.temp_obs)
    M=loadz(run+'/temp.npz')
    B=None
    if brun!=None and fexist(brun+'/temp.npz'): B=loadz(brun+'/temp.npz')
    plot_temp(M.bp.station,M,B,O,run=run,brun=brun,gd=gd,reftime=reftime,row=6,col=2,margin=[0.1,0.02,0.06,0.05],dxy=[0.08,0.00])

#%%
if p.harmonic: #only runable on cluster
    debug=True
    M=loadz(run+'/elevation.npz'); M.time=M.time+p.reftime
    B=None
    if brun!=None: B=loadz(brun+'/elevation.npz'); B.time=B.time+p.reftime
    O=loadz(p.elev_obs)
    E=loadz(p.elev_obs_info)
    lon={i:j for i,j in zip(E.station,E.lon)}
    lat={i:j for i,j in zip(E.station,E.lat)}
    station_name={i:j for i,j in zip(E.station,E.station_name)}
    stations=[i for i in M.bp.station if i.startswith('87')]
    if debug: stations=['8771341'] #['8779749','8775241','8773767','8772471','8771972','8771341','8770822','8764314','8760922','8735180']
    if not fexist(f'{run}/HA.npz') or 1:
        pstations,oHAs,mHAs,bHAs=[],[],[],[]
        for station in stations:#for harmonic analysis
            fp=(O.station==station)*(O.time>=M.time.min())*(O.time<=M.time.max())*(~isnan(O.wl))
            if sum(fp)==0: print('no obs data at',station); continue
            
            otimes=O.time[fp]
            owl=O.wl[fp]
            
            # get the largest continuous section
            ts = find_continuous_sections(otimes, 1.5/24)
            #print('largest continuous section',num2date(ts.msection),np.diff(ts.msection))
            fp=(otimes>=ts.msection[0])*(otimes<=ts.msection[1])
            otimes,owl=otimes[fp],owl[fp]
        
            #get model
            ist=(M.bp.station==station)
            fp=(M.time>=otimes.min())*(M.time<=otimes.max())
            mtimes,mwl=M.time[fp],M.elevation[ist,fp]
            
            if B!=None:
                bist=(B.bp.station==station)
                if sum(bist)==1:
                    fp=(B.time>=otimes.min())*(B.time<=otimes.max())
                    btimes,bwl=B.time,B.elevation[bist,fp].ravel()
                    bwl=bwl-mean(bwl)
            owl=owl-mean(owl)
            mwl=mwl-mean(mwl)
            #if len(mtimes)!=len(otimes): sys.exit('not consistent period between model and obs')
            oHA = harmonic_analysis(owl, median(diff(otimes)), otimes[0], tidal_names=p.tidal_names)
            mHA = harmonic_analysis(mwl, median(diff(mtimes)), mtimes[0], tidal_names=p.tidal_names)
            oHAs.append(oHA)
            mHAs.append(mHA)
            pstations.append(station)
            if B!=None: bHA = harmonic_analysis(bwl, median(diff(btimes)), btimes[0], tidal_names=p.tidal_names); bHAs.append(bHA)
            for tidal_name in ['K1','O1','Q1','P1','M2','S2','N2','K2']:
                itide=nonzero(p.tidal_names==tidal_name)[0][0]
                print(station,tidal_name,'observation vs model: {:.3f} vs {:.3f}, {:.0f} vs {:.0f}'.format(oHA.amplitude[itide],mHA.amplitude[itide],oHA.phase[itide]*180/pi,mHA.phase[itide]*180/pi))
            if debug:
                figure(figsize=[10,4])
                plot(otimes,owl,'k-')
                plot(btimes,bwl,'g-')
                plot(mtimes,mwl,'r-');
                legend(['model','obs']); show()
        S=zdata()
        S.oHAs,S.mHAs,S.bHAs,S.pstations=oHAs,mHAs,bHAs,pstations
        savez(f'{run}/HA.npz',S)
    else:
        S=loadz(f'{run}/HA.npz')
        oHAs,mHAs,bHAs,pstations=S.oHAs,S.mHAs,S.bHAs,S.pstations
    
    
    # plot the phase and amplitude at stations
    lons=array([lon[station] for station in pstations])
    lats=array([lat[station] for station in pstations])
    ind=in_domain(lons,lats,f'../{run}/hgrid.ll')
    for tidal_name in p.tidal_names:
        oa,ma,op,mp=[],[],[],[]
        itide=nonzero(p.tidal_names==tidal_name)[0][0]
        for oHA,mHA in zip(oHAs,mHAs):
            oa.append(oHA.amplitude[itide])
            ma.append(mHA.amplitude[itide])
            op.append(oHA.phase[itide])
            mp.append(mHA.phase[itide])
        oa,ma,op,mp=array(oa),array(ma),array(op),array(mp)
        p.gd.plot_bnd()
        
        subplot(211)
        p.gd.plot_bnd()
        fp=(ma>oa*1.1)*(ind==1)
        plot(lons[fp],lats[fp],'r.')
        fp=(ma<oa*0.9)*(ind==1)
        plot(lons[fp],lats[fp],'b.')
        title(run+': '+tidal_name)
        
        subplot(212)
        p.gd.plot_bnd()
        fp=((mp-op)%360<90)*(ind==1)
        plot(lons[fp],lats[fp],'r.')
        fp=((mp-op)%360>270)*(ind==1)
        plot(lons[fp],lats[fp],'b.')
        title(run+': '+tidal_name)
    #%%
    close('all')
    stations=['8779749','8775241','8773767','8772471','8771972','8771341','8770822','8764314','8760922','8735180']
    for m,tidal_name in enumerate(['K1','O1','M2']):
        itide=nonzero(p.tidal_names==tidal_name)[0][0]
        oa,ma,ba=[],[],[]
        for station in stations:
            fp=nonzero(array(S.pstations)==station)[0][0]
            oa.append(oHAs[fp].amplitude[itide])
            ma.append(mHAs[fp].amplitude[itide])
            ba.append(bHAs[fp].amplitude[itide])
        oa,ma,ba=array(oa),array(ma),array(ba)
        subplot(3,1,m+1)
        bar(arange(len(stations))-0.2,oa,width=0.2,color=[0.5,0.5,0.5,0.8])
        bar(arange(len(stations)),ba,width=0.2,color=[0.0,0.7,0.0,0.71])
        bar(arange(len(stations))+0.2,ma,width=0.2,color=[0.7,0.0,0.0,0.7])
        if m==0:legend(['Obs',brun,run])
        ylabel(tidal_name+' amp(m)')
    axes(position=[0.15,0.8,0.3,0.2])
    lons=[lon[station] for station in stations]
    lats=[lat[station] for station in stations]
    p.gd.plot_bnd()
    gca().axis('off')
    plot(lons,lats,'r.')
    
    for m,station in enumerate(stations):
        text(lons[m],lats[m],str(m),color='r')
    
    
    m2p_obs,m2p_mod,m2p_base=array(m2p_obs),array(m2p_mod),array(m2p_base)
    
    for itide,tidal_name in enumerate(p.tidal_names):
        amp_obs,amp_mod,amp_base,phase_obs,phase_mod,phase_base=[],[],[],[],[],[]
        for ist,station in enumerate(mstations):
            amp_obs.append(oHAs[ist].amplitude[itide+1])
            amp_mod.append(mHAs[ist].amplitude[itide+1])
            amp_base.append(bHAs[ist].amplitude[itide+1])
            phase_obs.append(oHAs[ist].phase[itide+1]*180/pi)
            phase_mod.append(mHAs[ist].phase[itide+1]*180/pi)
            phase_base.append(bHAs[ist].phase[itide+1]*180/pi)
        names=[station_name[i].split()[0] for i in mstations] #station name instead of NOAA ID
        phase_obs=array(phase_obs)
        phase_mod=array(phase_mod)
        phase_base=array(phase_base)

        #-- barplot for ampltiude and phase
        figure(figsize=[6,5])
        subplot(2,1,1)
        bar(arange(len(pstations))-0.2,amp_obs,0.2,color=[0,0,0,1])
        bar(arange(len(pstations)),amp_base,0.2,color=[0,1,1,0.5])
        bar(arange(len(pstations))+0.2,amp_mod,0.2,color=[1,0,0,0.5])
        legend(['Observation',brun,run,])
        ylabel(f'{tidal_name} amp. (m)')
        setp(gca(),xticks=arange(len(pstations)),xticklabels=names)
        xticks(rotation = 45)

        subplot(2,1,2)
        pdiff=(phase_base-phase_obs)%360; pdiff[pdiff>180]=pdiff[pdiff>180]-360
        bar(arange(len(pstations))-0.1,pdiff,0.2,color=[0,1,1,0.5])
        pdiff=(phase_mod-phase_obs)%360; pdiff[pdiff>180]=pdiff[pdiff>180]-360
        bar(arange(len(pstations))+0.1,pdiff,0.2,color=[1,0,0,0.5])
        ylabel(f'{tidal_name} phase diff (deg)')
        legend([f'{brun}-Obs',f'{run}-Obs'])
        setp(gca(),xticks=arange(len(pstations)),xticklabels=names)
        xticks(rotation = 45)

        tight_layout();
        figname=f'{run}/figs/wl_mainstem_harmonic_{tidal_name}.png'
        print(figname); savefig(figname,dpi=300)
