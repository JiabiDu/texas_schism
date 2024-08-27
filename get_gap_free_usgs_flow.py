#!/usr/bin/env python3
from pylib import *
close("all")
from pcom import *
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
# use stations downstream/upstream to fillup the gaps, use daily mean to find the best relation
#!! use local computer, since sklearn is not installed on sciclone
#%% Inuts
#some common inputs
begin_time=datenum(2000,1,1)
end_time=datenum(2024,6,1)
pretext,postext=[],[]
dll=2 #+- delta lon/lat
redo=False
if len(sys.argv)>1:  #in case this script is called from outside (e.g., for forecasting)
    begin_time=datenum(sys.argv[1])
    end_time=datenum(sys.argv[2])
    redo = sys.argv[3]=='True'

#% parameters for each river
rivers={'alabama':['02428400',['alabama','tombigbee','mobile','tensaw']],
        'tombigbee':['02469761',['alabama','tombigbee','mobile','tensaw']],
        'pearl':['02489500',['pearl','old pearl','mississippi']],
        'mississippi':['07374000',['mississippi','atchafalaya']], #use usgs
        'atchafalaya':['07381490',['mississippi','atchafalaya']], #temp since 2015
        'calcasieu':['08015500',['Calcasieu','sabine']],
        'sabine':['08030500',['sabine','neches','cow bayou']], #lacking 2012
        'neches':['08041780',['sabine','neches','cow bayou']],
        'west_fork_double_bayou':['08042558',['W Fk Double Bayou','trinity']], #; temp use Wallisville
        'trinity':['08067250',['trinity']], #near Romayor TX; temp use Wallisville
        'trinity2':['08067000',['trinity']], #trinity 2 at the liberty, TX
        'lost':['08067250',['trinity']], #no station, 10% of division is presumed (no ref, should find one)
        'cedar_bayou':['08067500',['cedar bayou','trinity']],
        'east_folk_goose_creek':['08067525',['goose','san jacinto']],
        'san_jancinto':['08072000',['san jancinto','lk houston']], #lake houston, water level will be used
        'greens_bayou':['08076000',['greens bayou','buffalo bayou','brays bayou']], 
        'buffalo_bayou':['08073600',['greens bayou','buffalo bayou','brays bayou']],
        'brays_bayou':['08075000',['greens bayou','buffalo bayou','brays bayou']],
        'sims_bayou':['08075400',['sims bayou','greens bayou','buffalo bayou','brays bayou']],
        #'clear_creek':['08077600',['clear','sims bayou','greens bayou','buffalo bayou','brays bayou']], #temp use Eagle Point
        'clear_creek2':['08076997',['clear','buffalo bayou']], #station at Pearland, less flow, less data gaps
        'chocolate_bayou':['08078000',['chocolate bayou','halls bayou']], 
        'brazos':['08116650',['brazos']], #temp since 2016
        'san_bernard':['08117705',['san bernard']],
        'colorado':['08162501',['colorado']],
        'tres_palacios':['08162600',['palacios','lavaca','carancahua']],
        'lavaca':['08164000',['lavaca','navidad','sandy']],
        'mission':['08189500',['mission','aransas']],
        'aransas':['08189700',['mission','aransas']],
        'guadalupe':['08188810',['guadalupe','san antonio']],
        'san_antonio':['08188060',['san antonio','guadalupe']], #merged into Guadalupe
        'nueces':['08211200',['nueces']], #into CC Bay
        'petronila_creek':['08212820',['petronila','san fernando','los olmos']], #into Baffin Bay
        'san_fernando_creek':['08211900',['petronila','san fernando','los olmos']], #into Baffin Bay
        'los_olmos_creek':['08212400',['petronila','san fernando','los olmos']], #into Baffin Bay
        'caloosahatchee_river':['02292900',['caloosahatchee',]], #tidal station
        'peace_river':['02296750',['peace',]],
        'myakka_river':['02298880',['myakka',]],
        'manatee_river':['02299950',['manatee',]], #including Gamble_creek (62cf
        'little_manatee_river':['02300500',['manatee','little manatee']], #temp using manatee
        'alafia_river':['02301500',['alafia',]],
        'hillsborough_river':['02304500',['hillsborough',]],
        'suwannee_river':['02323500',['suwannee',]],
        'steinhatchee_river':['02324000',['steinhatchee',]],
        'ochlockonee_river':['02330150',['ochlockonee',]],
        'sopchoppy_river':['02330150',['sopchoppy',]], #no measurement, using Ochlo
        'apalachicola_river':['02359170',['apalachicola',]], #big river
        'choctawhatchee_river':['02366500',['choctawhatchee',]],#big river
        'turkey_creek':['02250030',['turkey',]],
        'blackwater_river':['02370000',['blackwater',]],
        #'escambia_river':['02370500',['escambia',]],    #into Pensocola Bay
        'perdido_river':['02376500',['perdido',]],
        'fish_river':['02378500',['fish',]],
            }   
#%% parameters for each river
for iriver,river in enumerate(rivers):    
    print(iriver,'get river flow data for ',river)
    target_station=rivers[river][0]
    pretext=rivers[river][1]
    
    #%%
    y1=num2date(begin_time).year
    y2=num2date(end_time-1).year
    final_sname=f'npz/{target_station}_{y1}_{y2}' #gap free data
    if os.path.exists(final_sname+'.npz') and not redo: 
        print('exist '+final_sname+'.npz')
        continue
    S=loadz('stainfo.npz')
    fp=S.station==target_station
    print(S.staname[fp])
    
    # get stations
    station2={target_station:S.staname[fp]} # a dcit with staiton number and full name
    if len(pretext)>0:
        for ipre in pretext:
            tmp={i:j for i,j,lon,lat in zip(S.station,S.staname,S.lon,S.lat) if ipre.lower() in j.lower() and abs(lon-S.lon[fp])<dll and abs(lat-S.lat[fp])<dll}
            station2={**station2,**tmp}
    #print(station2) 
    stations=[i for i in station2]
    print(river,len(stations))
    
    #%% download the usgs data for the target station and stations within the same watershed and nearby stations
    sname=f'npz/usgs_{river}_{y1}_{y2}'
    if fexist(sname): print(f'{sname} exists already')
    if not fexist(sname) or redo: 
        get_usgs_flow(stations=stations,StartT=f'{y1}-1-1',EndT=f'{y2+1}-1-1',sname=sname,reRead=False,sdir=f'usgs_{y1}_{y2}',reDownload=False)
    
    #%% get daily data
    S=loadz(sname+'.npz')
    dailydata={}
    for station in unique(S.station):
        fp=S.station==target_station
        print('get daily mean for '+station)
        fp=S.station==station
        time1=S.time[fp]
        flow1=S.flow[fp]
        ltm=nanmean(flow1)
        time1,flow1=get_daily_mean(time1,flow1)
        dailydata[station]=[time1,flow1]
    #%% find best relation
    time1,flow1=dailydata[target_station]
    br2=[]
    b=[]
    for station in unique(S.station):
        time2,flow2=dailydata[station]
        best=linear_lag_analysis(time1,flow1,time2,flow2,daily=True,shiftings=arange(-5,5))
        br2.append(best.r2)
        b.append(best)
    br2=array(br2)
    b=array(b)
    bstation=unique(S.station)
    sind=argsort(br2)[::-1]
    br2=br2[sind]
    bstation=bstation[sind]
    b=b[sind]
    #%% plot the data
    [ps,pc]=get_subplot_position2(margin=[0.05,0.05,0.1,0.1],dxy=[0.05,0.00],ds=[6,2],dc=[0,0])
    ps=reshape(ps,(12,4)) #to make 3D dimension array to 2 dimension
    close('all')
    figure(figsize=[12,9])
    for m,station in enumerate(bstation):
        time1,flow1=dailydata[station]
        if m>11: break
        axes(position=ps[m])
        plot(time1,flow1,'ko',ms=1)
        xlim([datenum(2000,1,1),datenum(2022,1,1)])
        rtext(0.02,0.9,'R2={:.2f}; mean={:.2f}; {}={}'.format(br2[m],nanmean(flow1),station,station2[station]))
        grid('on')
    savefig(f'figs/{river}_original_data')
    
    [ps,pc]=get_subplot_position2(margin=[0.15,0.05,0.05,0.05],dxy=[0.05,0.05],ds=[3,3],dc=[0,0])
    ps=reshape(ps,(9,4)) #to make 3D dimension array to 2 dimension
    figure(figsize=[12,9])
    ip=-1
    for m,best in enumerate(b):
        pair1,pair2=best.pair1,best.pair2
        if len(best.pair1)==0: continue
        ip+=1
        if ip>8: break
        axes(position=ps[ip])
        plot(pair2,pair1,'k.',ms=1)
        plot(best.x_pred,best.y_pred,'r-')
        rtext(0.05,0.9,'R2={:.2f}'.format(best.r2))
        if ip==3: ylabel(target_station)
        rtext(0.1,0.3,bstation[m])
        rtext(0.1,0.1,station2[bstation[m]])
        grid('on')
    savefig(f'figs/{river}_linear_fits')
    #%% fitting with sklear linear regression
    # [ps,pc]=get_subplot_position2(margin=[0.15,0.05,0.05,0.05],dxy=[0.05,0.05],ds=[3,3],dc=[0,0])
    # ps=reshape(ps,(9,4)) #to make 3D dimension array to 2 dimension
    # ip=-1
    # figure(figsize=[9,8])
    # for m,best in enumerate(b):
    #     if best.R==1.0: continue
    #     ip+=1
    #     if ip>8: continue
    #     pair1,pair2=best.pair1,best.pair2
    #     x=pair2.reshape(-1,1)
    #     y=pair1
    #     x_pred=arange(min(x),max(x)).reshape(-1,1)
    #     transformer = PolynomialFeatures(degree=4, include_bias=False)
    #     transformer.fit(x)
    #     x_ = transformer.transform(x)
    #     model = LinearRegression().fit(x_, y)
    #     r_sq = model.score(x_, y)
    #     y_pred = model.predict(transformer.transform(x_pred))
        
    #     axes(position=ps[ip])
    #     plot(x,y,'k.',ms=1)
    #     plot(x_pred,y_pred,'r-')
    #     r_sq = model.score(x_, y)
    #     rtext(0.05,0.9,'R2={:.2f}'.format(r_sq))
    #     if ip==3: ylabel(target_station)
    #     rtext(0.3,0.1,bstation[m])
        
    
    #%% get gap-free data from interpolation
    #- for the missing values, interp from best related stations
    times=arange(begin_time+0.5,end_time) #when avg daily, 0.5 will be added
    time1=list(dailydata[target_station][0]) #convert to list to make it extentable
    flow1=list(dailydata[target_station][1]) #convert to list to make it extentable
    ltm=nanmean(flow1)
    #cut some of record to indicate the model's performance
    # time1=[]
    # flow1=[]
    for m,best in enumerate(b):
        missing_time=find_missing_time(time1, times)
        if m==0: omt=missing_time #original missing records
        if len(missing_time)==0: print('complete filling'); break
        if best.r2==1 or best.r2<0.8: continue
        print('Missing records before filling: {}'.format(len(missing_time)))    
        print('stations and R: ',bstation[m],best.r2)
        time2,flow2=dailydata[bstation[m]]
        time2=time2+best.shifting
        for imtime in missing_time:
            if imtime in time2:
                fp=time2==imtime
                if sum(fp)==1:
                    time1.append(imtime)
                    xp=best.transformer.transform(flow2[fp].reshape(-1,1))
                    yp=best.model.predict(xp)
                    flow1.append(yp[0])
    missing_time=find_missing_time(time1, times)
    print('Still missing records: {}'.format(len(missing_time)))
    time1=array(time1)
    flow1=array(flow1).astype('float')
    ind=argsort(time1)
    time1=time1[ind]
    flow1=flow1[ind]
    #%
    [xts,xls]=get_xtick(fmt=0)
    close('all')
    figure(figsize=[12,4])
    plot(time1,flow1,'-',color=[0,0,0,0.2])
    plot(dailydata[target_station][0],dailydata[target_station][1],'k.',ms=2)
    plot(omt,zeros(len(omt)),'bx',ms=1)
    if len(missing_time)!=0: plot(missing_time,ones(len(missing_time))*gca().get_ylim()[1],'o',color=[1,0,1,0.2],ms=1)
    setp(gca(),xticks=xts,xticklabels=xls,xlim=[begin_time,end_time])
    legend(['Final','Original','Mssing records','Missing records after interpolation'])
    ylabel('Flow (m3/s)')
    rtext(0.8,0.95,'Missing Records={}'.format(len(omt)))
    title('{}: {}'.format(target_station,station2[target_station]))
    tight_layout()
    savefig(f'figs/{river}_final_flow')
    
    if len(missing_time)>=60: 
        print('interp the gaps, number of missing records:{}'.format(len(missing_time)))
        flow1=interp(times,time1,flow1,left=ltm,right=ltm)
        time1=times
        plot(time1,flow1,'-',color=[0,0,0,0.5])
        savefig(f'figs/{river}_final_flow2')
    if len(missing_time)>0 and len(missing_time)<60: 
        print('interp scattered but still missing records')
        flow1=interp(times,time1,flow1); time1=times
    #if len(times)!=len(time1): sys.exit('missing records or there is replicate records')
    if len(time1)!=len(flow1) or len(time1)!=len(unique(time1)): sys.exit('there is replicate records')
    S=zdata()
    S.time=time1
    S.flow=flow1
    
    savez(final_sname,S)
    print('data saved into {}.npz'.format(final_sname))
