#!/usr/bin/env python3
from pylib import *
from pcom import *
#%%
def plot_land(gd,land_color=[0.7,0.7,0.7,1],ocean_color=[0,0.8,0.8,1]):
    xv,yv=gd.x[gd.ilbn[0]],gd.y[gd.ilbn[0]]
    xv=[*xv,min(xv),min(xv),max(xv)]
    yv=[*yv,min(yv),max(yv),max(yv)]
    fill(xv,yv,color=land_color)
    for ia,il in enumerate(gd.ilbn):
        if gd.island[ia]==1:fill(gd.x[il],gd.y[il],color=land_color)
    xlim([gd.x.min(),gd.x.max()])
    ylim([gd.y.min(),gd.y.max()])
    gca().set_facecolor(ocean_color)

def plot_salinity(stations,M,B,O,row,col,margin=[0.05,0.02,0.06,0.05],dxy=[0.03,0.00],p=None,ovar='seawater_salinity',ipt=1,figpre='',colw=6,
                  rowh=2,paper=False):
    run,brun,gd,reftime=p.run,p.brun,p.gd,p.reftime
    istack=0
    ps=get_subplot_position2(margin=margin,dxy=dxy,ds=[row,col]);ps=reshape(ps,(row*col,4))
    figure(figsize=[colw*col,2*rowh])
    xlims=[M.time.min()+reftime-2,M.time.max()*1.2+reftime]
    m=-1
    for station in stations:   
        fp=(O.station==station)*(O.var==ovar)*(O.time>M.time.min()+reftime)*(O.time<M.time.max()+reftime)
        if sum(fp)==0 and station!=stations[-1]: continue
        m+=1
        #plot the model surface and bottom
        print('plot salt for station ',m+1,station)
        if m%(col*row)==0:  clf()
        axes(position=ps[m%(col*row)])
        maxs,mins=0,35 #for ylim purpose
        
        olon,olat=O.lon[station],O.lat[station]
        fp=(O.station==station)*(O.var==ovar)*(O.time>=reftime)*(O.time<=M.time.max()+reftime)
        #plot(O.time[fp],O.data[fp],'.',color=[0.5,0.5,0.5,1],lw=1,ms=1)
        x,y=O.time[fp],O.data[fp]
        fp=~isnan(y)
        if sum(fp)>0:
            x,y=x[fp],y[fp]
            y=lpfilt(y, median(diff(x)), 0.48); y[-96:]=nan #get subtidal
            plot(x,y,'.',color=[0.1,0.1,0.1,1],lw=1,ms=1)
            maxs=max(maxs,max(y)); mins=min(mins,min(y))
            ox,oy=x,y
        
        #for bottom and surface in base run; only show the subtidal
        Bs=[B,]
        if type(B)==list: Bs=B
        for B in Bs: #multiple base runs
            if B!=None and sum(B.bp.station==station)>0:
                fp=(B.bp.station==station)*(B.bp.z==0); #surface
                fp=nonzero(fp)[0][0] 
                y=lpfilt(B.salinity[fp].ravel(), median(diff(B.time)), 0.48); y[-96:]=nan #get subtidal
                plot(B.time+reftime,y,'-',color=[0.3,0.3,1,1],lw=0.5)
                maxs=max(maxs,max(y)); mins=min(mins,min(y))
                Sb=get_skill(ox,oy,B.time+reftime,y)
                
                fp=(B.bp.station==station)*(B.bp.z==max(B.bp.z[B.bp.station==station])); #bottom
                fp=nonzero(fp)[0][0] 
                y=lpfilt(B.salinity[fp].ravel(), median(diff(B.time)), 0.48); y[-96:]=nan #get subtidal
                plot(B.time+reftime,y,'-',color=[0,0,0.7,1],lw=0.5)
                maxs=max(maxs,max(y)); mins=min(mins,min(y))

        # for surface
        fp=(M.bp.station==station)*(M.bp.z==0); fp=nonzero(fp)[0][0] #surface
        #plot(M.time+reftime,M.salinity[fp].ravel(),'-',color=[0,0.0,0.8,0.5],lw=0.2)
        # get the subtidal
        y=lpfilt(M.salinity[fp].ravel(), median(diff(M.time)), 0.48); y[-96:]=nan
        plot(M.time+reftime,y,'-',color=[1,0.3,0.3,1],lw=1)
        maxs=max(maxs,max(y)); mins=min(mins,min(y))
        Sm=get_skill(ox,oy,M.time+reftime,y)
        
        # for bottom
        fp=(M.bp.station==station)*(M.bp.z==max(M.bp.z[M.bp.station==station])); 
        fp=nonzero(fp)[0][0] #bottom
        #plot(M.time+reftime,M.salinity[fp].ravel(),'-',color=[0.8,0.0,0.0,0.5],lw=0.5)
        # get the subtidal
        y=lpfilt(M.salinity[fp].ravel(), median(diff(M.time)), 0.48); y[-96:]=nan
        plot(M.time+reftime,y,'-',color=[0.7,0.0,0,1],lw=1)
        maxs=max(maxs,max(y)); mins=min(mins,min(y))
        
        xlim(xlims) #times 1.2 to allow some space showing the monitoring station in a subset
        ylim([floor(mins),ceil(maxs)])
        set_xtick(fmt=0)
        rtext(0.02,0.9,'MAE:')
        rtext(0.1,0.9,'{:.0f} psu'.format(Sm.MAE),color='r')
        if B!=None: rtext(0.2,0.9,'{:.2f}'.format(Sb.MAE),color='b')
        
        if m%(row*col)<(row-1)*col and m!=len(stations)-1: setp(gca(),xticklabels=[]) #only the last row will have time stamp
        if m%(row*col)==0:ylabel('Salinity (psu)') 
        if not paper and m%(row*col)==ipt and B!=None: title(f'{run}(red); {brun}(blue); observation(gray)')
        if not paper and m%(row*col)==ipt and B==None: title(f'{run}(red); observation(gray)',fontsize=12,fontweight='bold')
        
        # show the monitoring station
        pa=gca().get_position()
        pa=[pa.x0,pa.y0,pa.width,pa.height]
        pa=[pa[0]+0.84*pa[2],pa[1]+pa[3]*0.05,pa[2]*0.15,pa[3]*0.9]
        axes(position=pa) 
        fill([olon-0.2,olon+0.2,olon+0.2,olon-0.2],[olat-0.4,olat-0.4,olat+0.4,olat+0.4],color=[1,1,1,1])
        #gd.plot_bnd(lw=1)
        #gO.plot(fmt=1,cmap='Blues',clim=[0,15],cb=False)
        #gd.plot_bnd(color='c')
        plot_land(gd)
        plot(olon,olat,'m.',ms=10)
        xlim([olon-0.2,olon+0.2]); ylim([olat-0.4,olat+0.4])
        rtext(0.1,0.9,station,fontweight='bold')
        
        setp(gca(),xticks=[],yticks=[])
        axis('off')
        gca().set_facecolor((1.0, 0.47, 0.42))
        if m%(row*col)==(row*col)-1 or station==stations[-1]:
            istack+=1
            #tight_layout()
            brunname=brun
            if type(brun)==list: brunname='_'.join(brun)
            figname=run+'/figs/salt{}_{:02d}_{}_vs_{}.png'.format(figpre,istack,run,brunname)
            print(figname); savefig(figname,dpi=300)

def get_skill(x1,y1,x2,y2):
    py1,py2=pair_data(x1,y1,x2,y2,hw=0.2/24)
    return get_stat(py1,py2,fmt=0)

#%%
def plot_salinity_one_layer(M,B,O,fmt=1,stations=['BOLI','MIDG','FISH','TRIN'],deps=[3,3,1.5,1.5],bay='Galveston',run=None,brun=None,gd=None,reftime=None,
                            runnames=None,titlename=None,mlg=0,show_location=False,row=None, col=None,margin=[0.1,0.02,0.06,0.05],dxy=[0.03,0.00],
                            ylims=None,dpi=300):
    lw=0.5; 
    if fmt==1: lw=1 #for subtidal signal
    Bs=[B,]
    if type(B)==list: Bs=B

    figure(figsize=[8,2*len(stations)]); offset=0;
    if row==None: row=len(stations); col=1
    ps=get_subplot_position2(margin=margin,dxy=dxy,ds=[row,col]);ps=reshape(ps,(row*col,4))

    ip=-1
    for m,station in enumerate(stations):
        fp=(O.station==station)*(O.var=='salinity')*(O.time>M.time.min()+reftime)*(O.time<M.time.max()+reftime)
        if sum(fp)==0: continue #no valid data is found
        ip+=1
        print('plot salt for station ',station) 

        axes(position=ps[ip])
        #plot observation
        ox,oy=O.time[fp],O.data[fp]
        if fmt==1: oy=lpfilt(O.data[fp],median(diff(O.time[fp])),0.48) #low pass filter
        plot(ox,oy,'x',color=[0,0,0,0.5],lw=lw,ms=1)
        
        #plot Model results from target run
        fp0=nonzero(M.bp.station==station)[0]
        tmpz=M.bp.z[fp0]
        tmp=abs(tmpz-deps[m]) # find the depth closest to the needed depth. 
        tfp=fp0[tmp==tmp.min()] # the station with closest depth 
        mx=M.time+reftime
        my=M.salinity[tfp].ravel()

        if fmt==1: my=lpfilt(M.salinity[tfp].ravel(),median(diff(M.time)),0.48) #0.48==50hour, lowppass filter
        plot(mx,my,'-',color=[1,0,0,0.7],lw=lw)    

        xlim([M.time.min()+reftime-2,M.time.max()+reftime])
        if show_location: xlim([M.time.min()+reftime-2,M.time.max()*1/0.85+reftime])
        if ylims==None: ylims=[0,36]
        ylim(ylims)
        set_xtick(fmt=2)  #show year ticks

        #get model skill
        S=get_skill(ox,oy,mx,my)

        bcolors=[[0,0.8,0.8,0.8],[0.8,0.,0.8,0.8],[0.0,0.8,0.0,0.8]]
        for k,B in enumerate(Bs): #multiple base runs
            if B!=None:
                # find the station with closest depth
                fp0=nonzero(B.bp.station==station)[0]
                tmpz=B.bp.z[fp0]
                tmp=abs(tmpz-deps[m])
                tfp=fp0[tmp==tmp.min()]
                y=B.salinity[tfp].ravel()
                if fmt==1: y=lpfilt(B.salinity[tfp].ravel(),median(diff(B.time)),0.48) #0.48==50hour
                plot(B.time+reftime,y,'-',color=bcolors[k],lw=lw*0.5)
                mxB=B.time+reftime
                myB=lpfilt(B.salinity[tfp].ravel(),median(diff(M.time)),0.48) #0.48==50hour

        if brun!=None: SB=get_skill(ox,oy,mxB,myB) #only show performance of the last base run
        if brun!=None and type(brun)!=list: 
            rtext(0.01,0.90,'('+chr(97+m)+') '+ r"$\bf{"+station+"}$"+'; MAE={:.2f} vs {:.2f} in {}'.format(S.MAE,SB.MAE,brun),backgroundcolor=[1,1,1,0.5])
        else:
            rtext(0.01,0.90,'('+chr(97+m)+') '+ r"$\bf{"+station+"}$"+'; MAE={:.2f} psu'.format(S.MAE),backgroundcolor=[1,1,1,0.5])
        grid('on',ls=':')

        if ip==mlg: #
            if not titlename==None: 
                title(titlename)
            if brun!= None: legend(['Obs.',run,brun])
            if brun==None: legend(['Obs.','Model'])
            if runnames!=None: legend(runnames)
        
        if ip!=len(stations)-1: setp(gca(),xticklabels=[])
        ylabel('Salinity')

        # show the monitoring station
        if show_location:
            olon,olat=O.lon[station],O.lat[station]
            pa=gca().get_position()
            pa=[pa.x0,pa.y0,pa.width,pa.height]
            pa=[pa[0]+0.85*pa[2],pa[1]+pa[3]*0.005,pa[2]*0.145,pa[3]*0.99]
            axes(position=pa) 
            #fill([olon-0.2,olon+0.2,olon+0.2,olon-0.2],[olat-0.4,olat-0.4,olat+0.4,olat+0.4],color=[1,1,1,1])
            plot_land(gd)
            gd.plot_bnd(c='k')
            plot(olon,olat,'m.',ms=10)
            xlim([olon-0.2,olon+0.2]); ylim([olat-0.4,olat+0.4])
            axis('off')

    #tight_layout()
    tmp=''
    if fmt==1: tmp='_sub_tidal'
    if type(brun)==list: brun='_'.join(brun)
    figname=f'{run}/figs/salt_{bay}{tmp}_{run}_vs_{brun}.png'
    print(figname); savefig(figname,dpi=dpi)
    
#%%
def plot_elevation(stations,M,B,O,row,col,margin=[0.1,0.02,0.06,0.1],dxy=[0.08,0.00],signal='full',
                   run=None,brun=None,reftime=None,station_name=None,xlims=None,ylims=None,gd=None,
                   regname='',paper=False,adj=True,lon=None,lat=None,show_location=True,dpi=300):
    #fmt=0: full signal
    #fmt=1: subtidal
    #fmt=2: tidal
    istack=0
    ps=get_subplot_position2(margin=margin,dxy=dxy,ds=[row,col]);ps=reshape(ps,(row*col,4))
    figure(figsize=[4.5*col,2*row])
    if xlims==None: xlims=[M.time.min()+reftime,M.time.max()+reftime]
    m=-1
    for station in stations:
        fp=(O.station==station)*(O.time>xlims[0])*(O.time<xlims[1])*(~isnan(O.wl))
        if sum(fp)==0 and station!=stations[-1]: continue
        m+=1
        if m%(col*row)==0:  clf()
        axes(position=ps[m%(col*row)])

        #new model data
        fp=M.bp.station==station
        foyi=lpfilt(M.elevation[fp].ravel(),diff(M.time).mean(),0.2)
        meanModel=mean(foyi)
        if xlims!=None and signal=='full': 
            meanModel=mean(M.elevation[fp].ravel()[(M.time+reftime<xlims[1])*(M.time+reftime>xlims[0])])
        mtime=M.time+reftime
        mfoyi=foyi
        mwl=M.elevation[fp].ravel() #full signal
        
        #old model
        if B!=None:   
            fp=B.bp.station==station
            bfoyi=lpfilt(B.elevation[fp].ravel(),diff(B.time).mean(),0.2)
            btime=B.time+reftime
            bwl=B.elevation[fp].ravel() #original full signal

        #observation
        fp=(O.station==station)*(O.time>xlims[0])*(O.time<xlims[1])*(~isnan(O.wl))
        offset=0
        if sum(fp)>0:
            otime=O.time[fp]
            owl=O.wl[fp]
            ofoyi=lpfilt(O.wl[fp],1/24,0.2)
            if adj==True: offset=nanmean(ofoyi)-meanModel
            #if xlims!=None and signal=='full':
            #    offset=nanmean(O.wl[fp])-meanModel #to correct the datum; TODO: use the median diff in subtidal signal
            print('offset=',offset)
            if signal=='full': plot(otime-reftime,owl,'-',color=[0,0,0,0.8],lw=0.5); 
            if signal=='subtidal': 
                plot(otime,ofoyi,'-',color=[0,0,0,1],lw=1)
                # get model performance
                Sm=get_skill(otime,ofoyi,mtime,mfoyi)
                if B!=None: Sb=get_skill(otime,ofoyi,btime,bfoyi)
            if signal=='tidal': #tidal signal 
                plot(otime-reftime,owl-ofoyi,'-',color=[0,0,0,0.6],lw=1); 
                if ylims==None:ylims=[min(owl-ofoyi)-0.2,max(owl-ofoyi)+0.2]
                ylim(ylims)
                
        #plot the base model results
        print('plot elevation for station ',m+1,station)
        if B!=None: #old model run
            fp=B.bp.station==station
            if signal=='full': plot(btime-reftime,bwl,'-',color=[0,0.8,0.8,0.6],lw=0.5) #full signal; use day instead date
            if signal=='subtidal': plot(btime,bfoyi,'-',color=[0.0,0.8,0.8,1],lw=1);btime,bfoyi=B.time+reftime,foyi
            if signal=='tidal': plot(btime-reftime,bwl-bfoyi,'-',color=[0.0,0.8,0.8,1],lw=0.5) #tidal signal
       
        #new model run
        if signal=='full': plot(mtime-reftime,mwl,'-',color=[1,0,0,0.6],lw=1) #full signal
        if signal=='subtidal': plot(mtime,mfoyi,'-',color=[0.8,0,0,1],lw=1);mtime,mfoyi=M.time+reftime,foyi #subtidal
        if signal=='tidal': plot(mtime-reftime,mwl-mfoyi,'-',color=[0.8,0,0,1],lw=1) #tidal signal

        # some annotations
        newxlims=xlims
        if show_location==True:
            newxlims=[xlims[0],xlims[1]+diff(xlims)[0]*1.2]
        if signal=='subtidal': xlim(newxlims); set_xtick(fmt=2)
        if signal=='tidal' or signal=='full':
           xlim(array(newxlims)-reftime)
           if m==col*(row-1): xlabel('Days since 2018-01-01 (GMT)')

        rtext(0.01,0.9,station+':'+station_name[station]) 
        if m%(col*row)==0: ylabel ('Water level (m)')
        if signal=='subtidal' and sum(fp)>0 and B!=None: rtext(0.01,0.8,'MAE: {:.2f}m vs {:.2f}m'.format(Sm.MAE,Sb.MAE))
        if signal=='subtidal' and sum(fp)>0 and B==None: rtext(0.01,0.8,'MAE: {:.2f}m'.format(Sm.MAE)) 
        if m%(row*col)<col*(row-1) and station != station[-1]: setp(gca(),xticklabels='')
        
        # show the monitoring station
        pa=gca().get_position()
        pa=[pa.x0,pa.y0,pa.width,pa.height]
        pa=[pa[0]+0.75*pa[2],pa[1]+pa[3]*0.01,pa[2]*0.25,pa[3]*0.98]
        axes(position=pa) 
        gd.plot_bnd(c='g',lw=1)
        olon,olat=lon[station],lat[station]
        plot(olon,olat,'mo',ms=10,markerfacecolor=[1,0,1,1])
        xlim([olon-0.2,olon+0.2]); ylim([olat-0.4,olat+0.4])
        fill([olon-0.2,olon+0.2,olon+0.2,olon-0.2],[olat-0.4,olat-0.4,olat+0.4,olat+0.4],color=[1,1,1,1])
        setp(gca(),xticks=[],yticks=[])
        axis('off')
        
        if m%(row*col)==(row*col)-1 or station==stations[-1]:
            istack+=1
            if signal=='full': figname=run+'/figs/elev_{}_{:02d}_{}_vs_{}.png'.format(regname,istack,run,brun)
            if signal=='subtidal': figname=run+'/figs/elev_subtidal_{}_{:02d}_{}_vs_{}.png'.format(regname,istack,run,brun)
            if signal=='tidal': figname=run+'/figs/elev_tidal_{}_{:02d}_{}_vs_{}.png'.format(regname,istack,run,brun)
            print(figname); savefig(figname,dpi=dpi)

#%% for current 
def get_direction(u,v):
    if v==0:
        if u>0: return 90
        if u<0: return 270
        if u==0: return 0
    else:
        direction=atan(u/v)*180/pi
        if v<0: direction=direction+180
    return direction%360

#% show the current station and name
def plot_obs_station(M,gd=None,var='current',p=None):
    run=p.run
    xlims=[M.bp.x.min()-0.2,M.bp.x.max()+0.2]
    ylims=[M.bp.y.min()-0.2,M.bp.y.max()+0.2]
    figure(figsize=[16,*8*diff(ylims)*2/diff(xlims)])
    if gd!=None: gd.plot(fmt=1,clim=[0,20],cb=False)
    plot(M.bp.x,M.bp.y,'ro')
    for lon,lat,station in zip(M.bp.x,M.bp.y,M.bp.station):
        text(lon,lat,station,fontsize=6)
    xlim(xlims)
    ylim(ylims)
    tight_layout()
    savefig(f'{run}/figs/{var}_obs_locations.png',dpi=300)

def plot_current_direction_hist(M,O,p=None):
    run=p.run
    figure(figsize=[20,9])
    for m,station in enumerate(M.bp.station): #M.bp.station:  #g06010 entrance, bin1 at the surface
        # plot the observational data
        ofp=O.station==station
        ubin=unique(O.bin[ofp])
        ofp=(O.station==station)*(O.bin==1)
        ox,oy=gappy_data(O.time[ofp],O.speed[ofp])
        subplot(5,7,m+1)
        hist(O.direction[ofp],bins=arange(5,360,10),color='k')
        mfp=M.bp.station==station
        direction=array([get_direction(u,v) for u,v in zip(M.horizontalVelX[mfp].ravel(),M.horizontalVelY[mfp].ravel())])
        hist(direction,bins=arange(5,360,10),color=[1,0.3,0.3,0.8])
        rtext(0.1,0.9,station+': {:.0f} vs {:.0f}'.format(median(O.direction[ofp]),median(direction)))
        #axis('off')
    tight_layout()
    savefig(f'{run}/figs/current_direction_hist.png',dpi=300)

def plot_current_speed_hist(M,O,p=None):
    run=p.run
    figure(figsize=[20,9])
    for m,station in enumerate(M.bp.station): #M.bp.station:  #g06010 entrance, bin1 at the surface
        # plot the observational data
        ofp=O.station==station
        ubin=unique(O.bin[ofp])
        ofp=(O.station==station)*(O.bin==1)
        subplot(5,7,m+1)
        hist(O.speed[ofp].ravel(),bins=arange(0,100,5),color='k')
        mfp=M.bp.station==station
        hist(M.speed[mfp].ravel()*100,bins=arange(0,100,5),color=[1,0.3,0.3,0.8])
        rtext(0.1,0.9,station+': {:.0f} vs {:.0f}'.format(median(O.speed[ofp]),median(M.speed[mfp]*100)))
        #axis('off')
    tight_layout()
    savefig(f'{run}/figs/current_speed_hist.png',dpi=300)

def plot_current(M,O,B,p=None,zoom=0,stations=None,row=5,col=1,reg=None,fmt=0,show_title=False,
                 run=None,brun=None,reftime=None,gd=None):
    #fmt: 0 for speed; 1 for east-west; 2 for south-north
    margin=[0.1,0.02,0.06,0.1]; dxy=[0.08,0.05]
    ps=get_subplot_position2(margin=margin,dxy=dxy,ds=[row,col]);ps=reshape(ps,(row*col,4))
    figure(figsize=[8*col,2*row])
    istack=0;m=-1
    if stations==None: stations=M.bp.station
    for station in stations: #M.bp.station:  #g06010 entrance, bin1 at the surface
        m+=1
        if m%(col*row)==0:  clf()
        axes(position=ps[m%(col*row)])
        
        # plot the observational data
        ofp=O.station==station
        ubin=unique(O.bin[ofp])
        colors=['k',[0.5,0.5,0.5,1]]
        #for n,ibin in enumerate([min(ubin),max(ubin)]): #the first and final bin
        for n,ibin in enumerate([min(ubin)+1,max(ubin)-1]): #the first and final bin
            ofp=(O.station==station)*(O.bin==ibin)
            if fmt==0: ox,oy=O.time[ofp],O.speed[ofp]
            if fmt==1: ox,oy=O.time[ofp],O.speed[ofp]*sin(pi/180*O.direction[ofp])
            if fmt==2: ox,oy=O.time[ofp],O.speed[ofp]*cos(pi/180*O.direction[ofp])
            ox,oy=gappy_data(ox,oy)
                
            if zoom!=0:  #zoom in the first few days (defined by zoom)
                fp=ox<=ox.min()+zoom
                ox,oy=ox[fp],oy[fp]
            plot(ox-reftime,oy,'-',color=colors[n],lw=0.5)
            
        #find the main direction, use histogram, or the medi
        if B!=None:
            mfp=B.bp.station==station
            tfp=B.time+reftime<=ox.max()
            if fmt==0: bx,by=B.time[tfp],B.speed[mfp,tfp].ravel()*100
            if fmt==1: bx,by=B.time[tfp],B.horizontalVelX[mfp,tfp].ravel()*100
            if fmt==2: bx,by=B.time[tfp],B.horizontalVelY[mfp,tfp].ravel()*100
            plot(bx,by,'-',color=[0.3,0.3,1,0.5],lw=1)
            
        mfp=M.bp.station==station
        tfp=M.time+reftime<=ox.max()
        if fmt==0: mx,my=M.time[tfp],M.speed[mfp,tfp].ravel()*100
        if fmt==1: mx,my=M.time[tfp],M.horizontalVelX[mfp,tfp].ravel()*100
        if fmt==2: mx,my=M.time[tfp],M.horizontalVelY[mfp,tfp].ravel()*100
        plot(mx,my,'-',color=[1,0.3,0.3,1],lw=1)
        xlim([ox.min()-reftime,ox.min()+(ox.max()-ox.min())*1.2-reftime])
        if fmt==0: ylim([0,150])
        if fmt==1 or fmt==2: ylim([-100,100])
        if station=='g06010' and fmt>0: ylim([-150,150])
        #set_xtick(fmt=1)
        rtext(0.05,0.9,station,backgroundcolor='w',fontweight='bold')
        
        #if m%(row*col)<(row-1)*col and m!=len(stations)-1: setp(gca(),xticklabels=[]) #only the last row will have time stamp
        velnames=['Current speed','Eastward speed','Northward speed']
        if m%(row*col)==0 and show_title: 
            title(f'{velnames[fmt]} {run}(red); {brun}(blue)')
            if B==None: title(f'{velnames[fmt]} {run}(red); {brun}(blue)')
            
        if m%(row*col)==(row-1)*col or station==M.bp.station[-1]: xlabel('Days since {}'.format(str(num2date(reftime))[0:10]))
        if m%(row*col)==0: 
            if brun==None: legend(['Obs. first bin','Obs. final bin','Model'],loc='upper center')
            if brun!=None: legend(['Obs. first bin','Obs. final bin',brun,run],loc='upper center')
            ylabel(velnames[fmt])
        
        #-- add station location
        # show the monitoring station
        olon,olat=M.bp.x[mfp],M.bp.y[mfp]
        pa=gca().get_position()
        pa=[pa.x0,pa.y0,pa.width,pa.height]
        pa=[pa[0]+0.84*pa[2],pa[1]+pa[3]*0.0,pa[2]*0.16,pa[3]*1.0]
        axes(position=pa) 
        gd.plot_bnd(lw=0.5)
        #gO.plot(fmt=1,cmap='Blues',clim=[0,15],cb=False)
        gd.plot_bnd(c='g')
        plot(olon,olat,'mo',ms=5,markerfacecolor=[1,0,1,1])
        xlim([olon-0.5,olon+0.5]); ylim([olat-0.5,olat+0.5])
        setp(gca(),xticks=[],yticks=[])
        axis('off')
        
        if m%(row*col)==(row*col)-1 or station==stations[-1]:
            istack+=1
            #tight_layout()
            xlabel('Days since {}'.format(str(num2date(reftime))[0:10]))
            figname=run+f'/figs/current_{reg}_fmt{fmt}_stack{istack}_zoom{zoom}_{run}_vs_{brun}.png'
            print(figname); savefig(figname)
#%%
def plot_shelf_current(O,M,B=None,run=None,brun=None,gd=None,stations=[42043,42050],station_nickname=['B','F'],xlims=None,
                       fmt=0,show_location=True,subtidal=False,reftime=datenum(2018,1,1),dpi=300):
    if xlims==None: xlims=[M.time.min(),M.time.max()]
    figure(figsize=[8,6])
    for ist,station in enumerate(stations):
        subplot(2,1,ist+1)
        # plot observation
        fp=(O.station==station)*(O.time>xlims[0]+reftime)*(O.time<xlims[1]+reftime)*(O.speed<=200)
        speed=O.speed[fp]; direction=O.direction[fp]; 
        if fmt==0: u=sin(direction*pi/180)*speed  #eastward speed
        if fmt==1: u=cos(direction*pi/180)*speed  #northward speed
        if subtidal: u=lpfilt(u,median(diff(O.time[fp])),0.48)
        plot(O.time[fp]-reftime,u/100,'k.',ms=2)

        if B!=None:
            fp=nonzero(B.bp.station==station_nickname[ist])[0][0]
            if fmt==0: x,y=B.time,B.horizontalVelX[fp] #eastward speed
            if fmt==1: x,y=B.time,B.horizontalVelY[fp] #northward speed
            if subtidal: y=lpfilt(y,median(diff(x)),0.48)
            x,y=x[x<xlims[1]],y[x<xlims[1]]
            plot(x,y,'-',lw=0.5,color=[0,1,1,0.5])

        # plot modeling
        fp=nonzero(M.bp.station==station_nickname[ist])[0][0]
        if fmt==0: x,y=M.time,M.horizontalVelX[fp] #eastward speed
        if fmt==1: x,y=M.time,M.horizontalVelY[fp] #northward speed
        if subtidal: y=lpfilt(y,median(diff(x)),0.48)
        x,y=x[x<xlims[1]],y[x<xlims[1]]
        plot(x,y,'-',lw=0.5,color=[1,0,0,0.5])

        # some annotations
        ylim([-1.5,1.5])
        if subtidal: ylim([-1,1])
        xlim([xlims[0],xlims[1]+diff(xlims)*0.2])
        if fmt==0: ylabel('Eastward speed (m/s)')
        if fmt==1: ylabel('Northward speed (m/s)')
        if subtidal==False: rtext(0.05,0.9,f'Buoy {station_nickname[ist]} - full')
        if subtidal: rtext(0.05,0.9,f'Buoy {station_nickname[ist]} - subtidal')
        if ist==1: xlabel('Days since 2018-01-01 (GMT)')

        # add monitoring location
        if show_location: 
            olon,olat=M.bp.x[fp],M.bp.y[fp]
            pa=gca().get_position()
            pa=[pa.x0,pa.y0,pa.width,pa.height]
            pa=[pa[0]+0.84*pa[2],pa[1]+pa[3]*0.0,pa[2]*0.16,pa[3]*1.0]
            axes(position=pa) 
            gd.plot_bnd(c='g')
            plot(olon,olat,'mo',ms=5,markerfacecolor=[1,0,1,1])
            xlim([olon-1,olon+1]); ylim([olat-0.2,olat+2.5])
            setp(gca(),xticks=[],yticks=[])
            axis('off')
    
    #save the figure
    tmp=''
    if subtidal: tmp='_subtidal'
    if fmt==0: tmp2='east'
    if fmt==1: tmp2='north'
    figname=run+f'/figs/shelf_current{tmp}_{tmp2}_{run}_vs_{brun}_d{int(round(xlims[0]))}_d{int(round(xlims[1]))}.png'
    print(figname); savefig(figname,dpi=dpi)
#%%
def plot_temp(stations,M,B,O,row,col,run=None,brun=None,gd=None,reftime=None,margin=[0.05,0.02,0.06,0.05],dxy=[0.03,0.00],ptitle=None):
    istack=0
    ps=get_subplot_position2(margin=margin,dxy=dxy,ds=[row,col]);ps=reshape(ps,(row*col,4))
    figure(figsize=[4*col,1.5*row])
    xlims=[M.time.min()+reftime-2,M.time.max()*1.2+reftime]
    m=-1
    for station in stations:   
        fp=(O.station==station)*(O.time>=reftime)*(O.time<=M.time.max()+reftime)
        if sum(fp)==0 and station!=stations[-1]: continue #no observation
        m+=1
        print('plot temperature for station ',m+1,station)
        if m%(col*row)==0:  clf()
        axes(position=ps[m%(col*row)])
        maxs,mins=0,35 #for ylim purpose
        olon,olat=M.bp.x[M.bp.station==station][0],M.bp.y[M.bp.station==station][0]
        
        x,y=O.time[fp],O.temp[fp]
        fp=~isnan(y)
        if sum(fp)>0:
            x,y=x[fp],y[fp]
            y=lpfilt(y, median(diff(x)), 0.48); y[-96:]=nan #get subtidal
            plot(x,y,'.',color=[0.5,0.5,0.5,1],lw=1,ms=1)
            maxs=max(maxs,max(y)); mins=min(mins,min(y))
        
        if B!=None and sum(B.bp.station==station)>0:
            fp=(B.bp.station==station)*(B.bp.z==0); #surface
            fp=nonzero(fp)[0][0] 
            y=lpfilt(B.temperature[fp].ravel(), median(diff(B.time)), 0.48); y[-96:]=nan #get subtidal
            plot(B.time+reftime,y,'-',color=[0.3,0.3,1,1],lw=0.5)
            maxs=max(maxs,max(y)); mins=min(mins,min(y))

        # for surface
        fp=(M.bp.station==station)*(M.bp.z==0); fp=nonzero(fp)[0][0] #surface
        #plot(M.time+reftime,M.salinity[fp].ravel(),'-',color=[0,0.0,0.8,0.5],lw=0.2)
        # get the subtidal
        y=lpfilt(M.temperature[fp].ravel(), median(diff(M.time)), 0.48); y[-96:]=nan
        plot(M.time+reftime,y,'-',color=[1,0.3,0.3,1],lw=1)
        maxs=max(maxs,max(y)); mins=min(mins,min(y))

        xlim(xlims) #times 1.2 to allow some space showing the monitoring station in a subset
        ylim([floor(mins),ceil(maxs)])
        set_xtick(fmt=0)
        if maxs<20:
            rtext(0.05,0.9,station,backgroundcolor='w')
        else:
            rtext(0.05,0.15,station,backgroundcolor='w')
        if m%(row*col)<(row-1)*col and m!=len(stations)-1: setp(gca(),xticklabels=[]) #only the last row will have time stamp
        if m%(row*col)==0: ylabel('Water temp (C)')
        if m%(row*col)==1 and B!=None and ptitle!=None: title(f'{run}(red); {brun}(blue); observation(gray)')
        if m%(row*col)==1 and B==None and ptitle!=None: title(f'{run}(red); observation(gray)',fontsize=12,fontweight='bold')
        
        # show the monitoring station
        pa=gca().get_position()
        pa=[pa.x0,pa.y0,pa.width,pa.height]
        pa=[pa[0]+0.85*pa[2],pa[1]+pa[3]*0.005,pa[2]*0.145,pa[3]*0.99]
        axes(position=pa) 
        fill([olon-0.2,olon+0.2,olon+0.2,olon-0.2],[olat-0.4,olat-0.4,olat+0.4,olat+0.4],color=[1,1,1,1])
        gd.plot_bnd(lw=1)
        #gO.plot(fmt=1,cmap='Blues',clim=[0,15],cb=False)
        gd.plot_bnd(color='c')
        plot(olon,olat,'m.',ms=10)
        xlim([olon-0.2,olon+0.2]); ylim([olat-0.4,olat+0.4])
        axis('off')
        
        if m%(row*col)==(row*col)-1 or station==stations[-1]:
            istack+=1
            #tight_layout()
            figname=run+'/figs/temp_{:02d}_{}_vs_{}.png'.format(istack,run,brun)
            print(figname); savefig(figname,dpi=300)

def remove_sticking(lon,lat):
    print('Removing sticking points')
    difflon=diff(lon,axis=0)
    difflat=diff(lat,axis=0)
    fp=(difflon==0)*(difflat==0)
    lon,lat=lon[1:,:],lat[1:,:]
    lon[fp]=nan; lat[fp]=nan
    return lon,lat

def remove_after_hitting_boundary(lon,lat,south=26,east=-87):
    print('Removing track after the particle hit the boundary')
    hit=(lat<=south)|(lon>=east)
    hits=sum(hit,axis=0)
    for ipar in arange(lon.shape[1]):
        if hits[ipar]==0: continue
        itime=nonzero(hit[:,ipar])[0][0] #the first hit
        lon[itime:,ipar]=nan
        lat[itime:,ipar]=nan
    return lon,lat

def cal_LET(lon,lat,xlims=[-100,-80],ylims=[18,31],xinterval=0.01,yinterval=0.01,time_interval=2,fname=None):
    print('Calculating LET'); t0=time.time()
    indi=np.round((lon-xlims[0])/xinterval).astype('int')
    indj=np.round((lat-ylims[0])/yinterval).astype('int')
    X,Y=meshgrid(arange(xlims[0],xlims[1],xinterval),arange(ylims[0],ylims[1],yinterval))
    Z=zeros(X.shape)
    for i,j,loni,lati in zip(indi.ravel(),indj.ravel(),lon.ravel(),lat.ravel()):
        if i<0 or j<0: continue #they are due to nan value
        Z[j,i]+=1 #be mindful it is [j,i],not [i,j]
    Z=Z/lon.shape[1]*time_interval  #Z/npar*time_interval
    print('Complete calculating LET in {:.2f}s'.format(time.time()-t0))
    if fname!=None: 
        S=zdata()
        S.X,S.Y,S.Z=X,Y,Z
        savez(fname,S); print(f'data saved into {fname}')
    return X,Y,Z

def get_age(lon,lat,time,rn=30*6): 
    #rn: number of release time for the simulation
    print('Calculating age')
    difflon=diff(lon,axis=0)
    difflat=diff(lat,axis=0)
    fp=(difflon!=0)*(difflat!=0)
    npar=lon.shape[1]
    age=zeros(lon.shape)
    rtime=zeros(npar)
    for ipar in arange(npar):
        itime=nonzero(fp[:,ipar])[0][0] #the first move
        rtime[ipar]=time[itime] #release time
    #use the median value for every npar/rn
    for itime in arange(rn):
        rtime[int(itime*npar/rn):int((itime+1)*npar/rn)]=median(rtime[int(itime*npar/rn):int((itime+1)*npar/rn)])
    for ipar in arange(npar):
        age[:,ipar]=time-rtime[ipar]
    age[age<0]=0
    return age,rtime

def remove_after_hitting_boundary(lon,lat,south=26,east=-87):
    print('Removing track after the particle hit the boundary')
    hit=(lat<=south)|(lon>=east)
    hits=sum(hit,axis=0)
    for ipar in arange(lon.shape[1]):
        if hits[ipar]==0: continue
        itime=nonzero(hit[:,ipar])[0][0] #the first hit
        lon[itime:,ipar]=nan
        lat[itime:,ipar]=nan
    return lon,lat

def cal_flushing(nc,reg=None,run=None,iplot=True,rn=30*6,south=26,east=-87):
    print('calculating flushing for',nc)
    import numpy as np
    C=ReadNC(nc)
    lon=ma.getdata(C.lon.val)
    lat=ma.getdata(C.lat.val)
    mtime=ma.getdata(C.time.val)/86400
    #calculate the age, find the release date
    t0=time.time()
    age,rtime=get_age(lon,lat,mtime,rn=rn)

    #remove those hiting southern boundary
    lon,lat=remove_after_hitting_boundary(lon,lat,south=south,east=east)

    print('Cal inpolygon') #take about 20s
    bp=read_schism_reg(reg)
    ind_inpoly=reshape(inside_polygon(c_[lon.ravel(),lat.ravel()],bp.x,bp.y)==1,lon.shape)

    # get release index
    print('Get release index')
    ind_release=tile(np.round((rtime-rtime[0])/mean(diff(unique(rtime)))).astype('int'),[lon.shape[0],1])
    if rn==1: ind_release=tile(0,[lon.shape[0],1])

    # get age index
    print('Get age index')
    ind_age=np.round((age/median(diff(age[:,0])))).astype('int')

    # get the number inside the polygon
    print('Get number of particle in the polygon over time for each release'); 
    Z=zeros([rn,len(unique(age[:,0]))])
    for i,j,k in zip(ind_release.ravel(),ind_age.ravel(),ind_inpoly.ravel()):
        if i<0 or j<0: continue
        if k: Z[i,j]+=1
    print('Complete calculating in {:.2f}s'.format(time.time()-t0))
    lage=unique(age[:,0])
    npar=lon.shape[1]/rn
    S=zdata()
    S.lage=lage; S.release_time=unique(rtime)
    S.Z=Z
    S.ind_release=ind_release
    S.ind_age=ind_age
    S.ind_inpoly=ind_inpoly
    savez('tmp.npz',S)

    from sklearn.linear_model import LinearRegression
    lm = LinearRegression()
    flushing=[]
    efolding=[]
    print('Caculating flushing time')
    for ir in arange(rn):
        x=lage[1:90*12].reshape(-1, 1)
        y=Z[ir,1:90*12]/npar
        x=x[y>0]; y=y[y>0]
        lm.fit(x,log(y))
        flushing.append(-1/lm.coef_)
        if sum(Z[ir,:]/npar<exp(-1))>0: 
            efolding.append(lage[nonzero(Z[ir,:]/npar<exp(-1))[0][0]])
        else:
            efolding.append(nan)

    S=zdata()
    S.lage=lage; S.release_time=unique(rtime)
    S.Z=Z
    S.npar=npar
    S.flushing=array(flushing)
    S.efolding=array(efolding)
    savez('npz/'+run+'_flushing.npz',S)
    
    if iplot: 
        close('all')
        print('Plotting')
        figure(figsize=[6,4])
        for ir in arange(rn):
            plot(lage[1:],Z[ir,1:]/npar*100,lw=0.5,color=[0,0,0,0.2])
        plot([0,90],[exp(-1)*100,exp(-1)*100],'k--')
        plot(lage[1:],Z[:,1:].mean(axis=0)/npar*100,lw=1.5,color=[1,0,0,1])
        xlim([0,90])
        ylabel('Percent of paticles retained (%)')
        xlabel('Days since release')
        rtext(0.2,0.83,'Mean flushing time (from LR)= {:.2f} \u00B1 {:.2f} d'.format(mean(flushing),std(flushing)))
        rtext(0.2,0.9,'Mean flushing time (from e-folding)= {:.2f} \u00B1 {:.2f} d'.format(mean(efolding),std(efolding)))
        tight_layout()
        grid('on')
        savefig('figs/'+run+'_flushing.png',dpi=300)


