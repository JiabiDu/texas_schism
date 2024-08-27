#!/usr/bin/env python3
from pylib import *
close("all")
from pcom import *
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
#%% 
begin_time=datenum(2023,1,1)
end_time=datenum(2024,5,17)
reRead=False
npz_dir='npz_2024'
os.makedirs(npz_dir,exist_ok=True)

#%% get water temperature data 
y1,y2=num2date(begin_time).year,num2date(end_time-1).year
         #miss        atcha      trinity      neuces
stations=['07374000','07381600','08067252', '08211200','02292900','02299230','023000095','023000095','02301718','02306028','02323592','02323592','02327031','02359170','02327031','02359170',]
get_usgs_temp(stations=stations,StartT=f'{y1}-1-1',EndT=f'{y2+1}-1-1',sname=f'{npz_dir}/usgs_temp_{y1}_{y2}',reRead=reRead,sdir=None)
U=loadz(f'{npz_dir}/usgs_temp_{y1}_{y2}.npz')
figure(figsize=[16,10])
for m,station in enumerate(stations):
    fp=U.station==station
    subplot(4,4,m+1)
    plot(U.time[fp],U.temp[fp])
    xlim([datenum(y1,1,1),datenum(y2+1,1,1)])
    set_xtick()
tight_layout()
savefig('figs/usgs_temp')

         #albama        pearl          calcieu        sabine         san Jacinto     clear ck       chocolate bayou colorado       lavaca       mission        gurdalupe     baffin
stations=['NOAA8737048','NOAA8747437', 'NOAA8767961', 'NOAA8770475', 'NOAA8770777', 'NOAA8771013', 'NOAA8771972', 'NOAA8773146', 'NOAA8773259','NOAA8774770', 'NOAA8773037', 'NOAA8776604']
stations=[i.replace('NOAA','') for i in stations]
get_noaa_tide_current(stations=stations,years=arange(y1,y2+1),varnames=['water_temperature'],sdir=f'noaa_{y1}_{y2}/')
process_noaa_tide_current(stations=stations,years=arange(y1,y2+1),varnames=['water_temperature'],sname_pre=f'{npz_dir}/noaa_',sdir=f'noaa_{y1}_{y2}/',read_again=reRead)

N=loadz(f'{npz_dir}/noaa_water_temperature_{y1}_{y2}.npz')
figure(figsize=[12,8])
for m,station in enumerate(stations):
    fp=N.station==station
    subplot(6,2,m+1)
    plot(N.time[fp],N.temp[fp])
    xlim([datenum(2000,1,1),datenum(2022,1,1)])
    set_xtick()
tight_layout()
savefig('figs/noaa_temp')

S=zdata()
fp=isnan(N.temp)
S.station=concatenate([U.station,N.station[~fp]])
S.time=concatenate([U.time,N.time[~fp]])
S.temp=concatenate([U.temp,N.temp[~fp]])
#%% get daily data
close('all')
figure(figsize=[15,10])
dailydata={}
for m,station in enumerate(unique(S.station)):
    subplot(6,4,m+1)
    print('get daily mean for '+station)
    fp=S.station==station
    time1=S.time[fp]
    flow1=S.temp[fp]
    ltm=nanmean(flow1)
    time1,flow1=get_daily_mean(time1,flow1,lpf=False)
    plot(time1,flow1,'r-',lw=0.5)
    time1,flow1=remove_spike(time1,flow1,ds=2,hw=5,rmnan=True) #remove spiking value and remove nan values
    plot(time1,flow1,'k-',lw=0.5)
    set_xtick()
    
    #get seasonal mean
    print('get monthly mean')
    mon=[num2date(i).month for i in time1]
    mmean,mstd,mmin,mmax,mtime=arange(12)*1.0,arange(12)*1.0,arange(12)*1.0,arange(12)*1.0,arange(12)*1.0
    for imon in arange(12):
        tmp=flow1[mon==imon+1]
        mmean[imon]=tmp.mean()
        mstd[imon]=tmp.std()
        mmin[imon]=tmp.min()
        mmax[imon]=tmp.max()
        mtime[imon]=datenum(2000,imon+1,15)
        
    dailydata[station]=[time1,flow1,mmean,mstd,mmin,mmax] 
    if sum(isnan(flow1))>0: sys.exit('nan value exists')
        
figure(figsize=[15,10])
for m,station in enumerate(unique(S.station)):
    time1,flow1,mmean,mstd,mmin,mmax=dailydata[station]
    subplot(6,4,m+1)
    plot(mtime,mmean)
    plot(mtime,mmean+mstd)
    plot(mtime,mmean-mstd)
    plot(mtime,mmin,'--')
    plot(mtime,mmax,'--')
    title(station)
    if m==0:
        legend(['mean','+std','-std','min','max'])
tight_layout()
savefig('figs/river_seasonal')
#%% for each station find the full record
for target_station in unique(S.station):
    close('all')
    final_sname=f'{npz_dir}/{target_station}_{y1}_{y2}' #gap free data
    
    #% find best relation and rank them with r2
    time1,flow1,mmean,mstd,mmin,mmax=dailydata[target_station]
    br2=[]
    b=[]
    for station in unique(S.station):
        time2,flow2,*tmp=dailydata[station]
        best=linear_lag_analysis(time1,flow1,time2,flow2,daily=True,shiftings=[0],degree=2) #flow1 is the target
        br2.append(best.r2)
        b.append(best)
    br2=array(br2)
    b=array(b)
    bstation=unique(S.station)
    sind=argsort(br2)[::-1]
    br2=br2[sind]
    bstation=bstation[sind]
    b=b[sind]
    
    #% plot the data
    [ps,pc]=get_subplot_position2(margin=[0.15,0.05,0.05,0.05],dxy=[0.05,0.05],ds=[3,3],dc=[0,0])
    ps=reshape(ps,(9,4)) #to make 3D dimension array to 2 dimension
    figure(figsize=[12,9])
    ip=-1
    for m,best in enumerate(b):
        pair1,pair2=best.pair1,best.pair2
        ip+=1
        if ip>8: break
        axes(position=ps[ip])
        plot(pair2,pair1,'k.',ms=1)
        plot(best.x_pred,best.y_pred,'r-')
        rtext(0.05,0.9,'R2={:.2f}'.format(best.r2))
        if ip==3: ylabel(target_station)
        rtext(0.1,0.3,bstation[m])
        grid('on')
    savefig(f'figs/{target_station}_linear_fits')
    
    #% get gap-free data from interpolation
    #- for the missing values, interp from best related stations
    times=arange(begin_time+0.5,end_time) #when avg daily, 0.5 will be added
    time1=list(dailydata[target_station][0]) #convert to list to make it extentable
    flow1=list(dailydata[target_station][1]) #convert to list to make it extentable
    ltm=nanmean(flow1)
    #cut some of record to indicate the model's performance
    # time1=[]
    # flow1=[]
    for m,best in enumerate(b):
        missing_time=find_missing_time(time1, times,block=5)
        if m==0: omt=missing_time #original missing records
        if len(missing_time)==0: print('complete filling'); break
        if best.r2==1 or best.r2<0.85: continue
        mdff=(best.pair1-best.pair2).mean()
        print('Missing records before filling: {}'.format(len(missing_time)))    
        print('stations, r2, mean diff: ',bstation[m],best.r2,mdff)
        time2,flow2,*tmp=dailydata[bstation[m]]
        time2=time2+best.shifting
        
        #% simply use the observation record, but add the mean difference 
        # for imtime in missing_time:
        #     if imtime in time2:
        #         fp=time2==imtime
        #         if sum(fp)==1:
        #             time1.append(imtime)
        #             flow1.append(flow2[fp]+mdff)
        #%use the linear relationship
        for imtime in missing_time:
            if imtime in time2:
                fp=time2==imtime
                if sum(fp)==1:
                    time1.append(imtime)
                    xp=best.transformer.transform(flow2[fp].reshape(-1,1))
                    yp=best.model.predict(xp)
                    flow1.append(yp[0])
    missing_time=find_missing_time(time1, times,block=3)
    print('Still missing records: {}'.format(len(missing_time)))
    time1=array(time1)
    flow1=array(flow1).astype('float')
    ind=argsort(time1)
    time1=time1[ind]
    flow1=flow1[ind]
    
        
    if len(missing_time)>=60: 
        print('interp the gaps, number of missing records:{}'.format(len(missing_time)))
        flow1=interp(times,time1,flow1,left=ltm,right=ltm)
        time1=times
    if len(missing_time)>0 and len(missing_time)<60: 
        print('interp scattered but still missing records')
        flow1=interp(times,time1,flow1); time1=times
    if len(times)!=len(time1): sys.exit('missing records or there is replicate records')
    
    #remove spiking values
    time1,flow1=remove_spike(time1,flow1,ds=2,hw=5,inter_spike=True) #remove spiking value and remove nan values
    
    #based on observed seasonality, thoese exceeding the min and max will be replaced with min and max
    time2,flow2,mlower,mupper=seasonal_impose(time1,flow1,mtime,mmean,mstd,mmin,mmax,ds=None)
    
    #% plot the final result
    [xts,xls]=get_xtick(fmt=0)
    figure(figsize=[12,4])
    plot(time1,flow1,'-',color=[0,0,0,0.2],lw=1)
    plot(time2,flow2,'r-',lw=0.5)
    plot(time2,mlower,'--',lw=0.1)
    plot(time2,mupper,'--',lw=0.1)
    plot(dailydata[target_station][0],dailydata[target_station][1],'k.',ms=2)
    plot(omt,zeros(len(omt)),'bx',ms=2)
    plot(missing_time,ones(len(missing_time))*gca().get_ylim()[1],'o',color=[1,0,1,0.2],ms=2)
    setp(gca(),xticks=xts,xticklabels=xls,xlim=[begin_time,end_time])
    legend(['Final','Final2','lower','upper','Original','Mssing records','Missing records after interpolation'])
    ylabel('Temp (C)')
    rtext(0.8,0.95,'Missing Records={}'.format(len(omt)))
    title('{}'.format(target_station))
    tight_layout()
    savefig(f'figs/final_temp_{target_station}')
    
    Z=zdata()
    Z.time=time2
    Z.temp=flow2
    
    print('final data saved into '+final_sname+'.npz')
    savez(final_sname,Z)
    
