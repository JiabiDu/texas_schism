#!/usr/bin/env python3
import dash
# import dash_html_components as html
# import dash_core_components as dcc
from dash import html,dcc
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import pandas as pd
from datetime import datetime
import numpy as np
import plotly.express as px
from pylib import *
from scipy.signal import butter, filtfilt

#% 
default_station='8771450'
run='run19p'
brun=None
ranges=[[20,50],[0,730],[20,50]]
adj=False

#%
M=loadz(f'{run}/elevation.npz')
if brun!=None: B=loadz(f'{brun}/elevation.npz')
O=loadz('hourly_height_2000_2022.npz'); O.time=O.time-datenum(2018,1,1)
E=loadz('stainfo.npz')
lon={i:j for i,j in zip(E.station,E.lon)}
lat={i:j for i,j in zip(E.station,E.lat)}
station_name={i:j for i,j in zip(E.station,E.station_name)}

#% 
def lpf(data,normal_cutoff=1/24): #low pass filter
    order=2
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

def gappy_data(x,y,gfactor=10):
    x,y=array(x),array(y)
    if len(x)==0: return array([])
    ind=argsort(x)
    x,y=x[ind],y[ind]
    mdt=median(diff(x))
    tmpx=[]
    for m,xi in enumerate(x):
        if m==0: continue
        dt=xi-x[m-1]
        if dt>mdt*gfactor: #gap larger than multiple times of the median dt^M
           tmpx.append(x[m-1]+mdt)
    x=concatenate([x,array(tmpx)])
    y=concatenate([y,ones(len(tmpx))*nan])
    ind=argsort(x)
    x,y=x[ind],y[ind]
    return x,y
def plot_water_level(station,fmt=0,adj=False):
    fp=nonzero(M.bp.station==station)[0][0]
    y=M.elevation[fp].ravel(); mwl=y
    if fmt>=1: wl_sub=lpf(y,1/25)    
    if fmt==1: y=wl_sub #subtidal
    if fmt==2: y=y-wl_sub #tidal
    fig_ts=px.line(x=M.time,y=y)
    fig_ts.update_traces(line=dict(color='red', width=2),name=run,showlegend=True)
    
    if brun!=None:
        fp=nonzero(B.bp.station==station)[0][0]
        y=B.elevation[fp].ravel(); 
        if fmt>=1: wl_sub=lpf(y,1/25)    
        if fmt==1: y=wl_sub #subtidal
        if fmt==2: y=y-wl_sub #tidal
        fig_ts=fig_ts.add_scatter(x=B.time,y=y,mode='lines',showlegend=True,line=dict(color='green', width=1),name=brun)
    
    # add observational data
    fp=(O.station==station)*(O.time>0)*(O.time<M.time.max())
    if sum(fp)==0: return fig_ts
    otime,owl=gappy_data(O.time[fp],O.wl[fp])
    offset=nanmean(mwl)-nanmean(owl)
    if adj: owl+=offset
    y=owl
    ttitle=f'Water level at {station} {station_name[station]}'; xtitle='Days since 2018-1-1'
    if fmt>=1: wl_sub=lpf(owl,1/50);  xtitle=''
    if fmt==1: y=wl_sub; ttitle='Subtidal water level'; #subtidal
    if fmt==2: y=y-wl_sub; ttitle='Tidal water level'; #tidal

    fig_ts.add_scatter(x=otime,y=y,mode='lines',showlegend=True,line=dict(color='black', width=1),name='Observation')
    fig_ts.update_layout(yaxis=dict(range=[-1.,1.],title='Water level (m)'),
                          xaxis=dict(range=ranges[fmt],title=xtitle),
                          title=ttitle,title_x=0.5)
    return fig_ts



# add the map
fig=px.scatter_geo(lon=M.bp.x,lat=M.bp.y,hover_name=M.bp.station)
fig = px.scatter_mapbox(lon=M.bp.x,lat=M.bp.y,hover_name=M.bp.station, zoom=3)
fig.update_layout(clickmode='event+select', width=500,height=800)
fig.update_traces(marker_size=10,marker_color='darkblue')
fig.update_layout(mapbox_style="open-street-map", mapbox_zoom=4, mapbox_center_lat = 29, mapbox_center_lon=-94)

# get data
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

#% dashboard layout
app = dash.Dash(external_stylesheets =external_stylesheets)   #initialising dash app
app.layout = html.Div([
    html.H1("TxSCHISM Water Level comparison",style={'textAlign': 'center'}),  # Header
    # html.Div([
    #     html.Label("run:  ", style={'fontSize': '30px'}),  # Label for color input
    #     dcc.Input(id='run', type='text', value='run19p', style={'fontSize': '20px', 'width':'150px','height':'30px'}),  # Input box for color selection
    #     html.Label("base run:  ", style={'fontSize': '30px'}),  # Label for data variable dropdown
    #     dcc.Input(id='brun', type='text', value=brun, style={'fontSize': '20px', 'width':'150px','height':'30px'}),  # Dropdown menu for data variable selection
    # ], style={'textAlign': 'center'}),  # Flex display to arrange items in a row
    html.Div([
    dcc.Graph(id='map',figure=fig),
    html.Div([
            dcc.Graph(id='waterlevel_plot', figure=plot_water_level(default_station,adj=adj)),  # Water level plot on the top right
            dcc.Graph(id='subtidal_plot', figure=plot_water_level(default_station,fmt=1,adj=adj)),  # Other plot on the bottom right
            dcc.Graph(id='tidal_plot', figure=plot_water_level(default_station,fmt=2,adj=adj)),
        ], style={'display': 'flex', 'flexDirection': 'column', 'flex': '1'})  # Flexbox for right column
    ], style={'display': 'flex'})
])

@app.callback(
    [Output('waterlevel_plot', 'figure'),
     Output('subtidal_plot', 'figure'),
     Output('tidal_plot', 'figure')],
    [Input('map', 'clickData')]
)
def update_timeseries(clickData):
    if clickData is not None:
        station = clickData['points'][0]['hovertext']
        return plot_water_level(station,adj=adj),plot_water_level(station,fmt=1,adj=adj),plot_water_level(station,fmt=2,adj=adj)
    else:
        return plot_water_level(default_station)
    
if __name__ == '__main__':  #after run the entire scrip; go to http://127.0.0.1:8050
    app.run_server()