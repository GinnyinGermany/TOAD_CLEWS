# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 12:35:38 2021

@author: impy2
"""
from pathlib import Path
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
# plt.rcParams.update({'font.size': 24})
from scipy.io import loadmat

################# Specify the experiment and models ################

experiment = '1pctCO2'
variant_id = 'r1i1p1f1'

################# Specify variable of interest ######################
# var2 is the same as var but all lower case
var = 'cVeg' # or treeFrac
var2 = 'cveg' # or treefrac

############## Specify units of variable ###########################
units = '$kgm^{-2}$'

#create subplots
fig, axes = plt.subplots(1, 1, figsize=(6, 5.5), sharex=True)

#define models to analyse and relevant latitudes and longitudes for each
models = ['EC-Earth3-Veg', 'GFDL-ESM4', 'MPI-ESM1-2-LR', 'NorCPM1', 'SAM0-UNICON', 'TaiESM1']
lats = [-1,-5, -14, 0, 1, 0]
lons = [-69, -65, -72, -60, -55, -60]
colors = ['tab:blue', 'tab:green', 'tab:purple', 'darkred', 'dimgray', 'goldenrod']

#loop through each model
for i in range (0,6):

    model = models[i]
    latitude = lats[i]
    longitude = lons[i]

    if model == 'EC-Earth3-Veg':
    
        date_range = '1850-2000'
        
    elif model == 'GFDL-ESM4':
        date_range = '1850-1999'
        
    elif model == 'MPI-ESM1-2-LR':
        
        date_range = '1850-2014'
    
    elif model == 'NorCPM1':
        
        date_range = '1850-2013'
    
    elif model == 'TaiESM1':
        
        date_range = '1850-2000'
    
    elif model == 'SAM0-UNICON':
        
        date_range = '1850-1999'
        
    elif model == 'UKESM1-0-LL':
        date_range = '1850-2000'
        
    ##################### Set RELDIFF True if taking relative differennce #################
    ##################### False if absolute difference ####################################    
    RELDIFF = False
    
    # Is variable measured in per second?
    PERSEC = False
        
    ############ Specify name of directory of data which you want to plot #################
    region = 'Amazon'
    
    ############ Specify name of directory of that stores the figures #################
    region2 = 'Amazon'
    
    # Window length used to calculate abrupt shift over
    wl = 15
    
    ################## Specify path to processed data directory #####################
    path = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment+'/Processed_data_monthly/'+region+'/'
    path2 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment+'/Analysis_data/'+region2+'/'

    # File name of interpolated data
    fname = path+var+'_'+model+'_'+experiment+'_'+variant_id+'_'+date_range+'_mon_'+region+'.nc'
    fname = Path(fname)
    
    # # Load in data
    f = nc.Dataset(fname,'r')
    mat = loadmat(path2+var+'_'+model+'_as_grads_wl'+str(wl)+'_data_'+experiment+'_'+region2+'.mat')
    
    # #extract indicies of lon and lat in the matrices 
    lon = np.array(f.variables['lon'])
    lon_point = np.where(lon==longitude)
    
    lat = np.array(f.variables['lat'])
    lat_point = np.where(lat==latitude)
    
    # Extract data and if measured in per second convert to per year
    if PERSEC:
        x = f.variables[var2][:]*3600*24*360
    else:
        x = f.variables[var2][:]
    
    # Close dataset
    f.close()
      
    ### extract the data as a time series ####
    data = x[:, lat_point, lon_point]
    data = data.flatten()
    #get rid of last datapoint for NorCPM1
    data=data[0:-1]
    t = np.arange(1, data.size+1)
    t2 = t/12
    
    # Create array of CO2 for 1pctCO2 run
    co2_start = 284.3186666723341
    co2 = np.zeros(data.size)
    for k in range(data.size):
        co2[k] = co2_start*(1.01**k)
    
    ##extracting the CO2 conc and time point at which abrupt shift is detected 
    co2_tip = mat['co2_tip']
    co2_at_tip = co2_tip[lat_point, lon_point]
    tip_index = np.where(co2==co2_at_tip)
    time_at_tip = t[tip_index[1]]                     

    #define functions for creating secondary x axis
    def forward(x):
        return co2_start*(1.01**x)
    
    # co2 to time
    def inverse(x):
        return np.log(x/co2_start)/np.log(1.01)
    
    ############## plot the time series data ############
    if latitude>0:
        lat_str = str(int(np.linalg.norm(latitude)))+u'\N{DEGREE SIGN}N'
    else:
        lat_str = str(int(np.linalg.norm(latitude)))+u'\N{DEGREE SIGN}S'
        
    lon_str = str(int(np.linalg.norm(longitude)))+u'\N{DEGREE SIGN}W'
    
    if model == 'NorCPM1' or model == 'MPI-ESM1-2-LR':
        axes.plot(t2[0:1799], data[0:1799], label=model, color=colors[i])
    else:
        axes.plot(t2, data, label=model, color=colors[i])
    
plt.xlabel('Time (model years)')
plt.ylabel('cVeg ('+units+')')
plt.xticks(np.arange(0, 200, 50))
plt.yticks(np.arange(0, 40, 10))
secax = axes.secondary_xaxis('top', functions=(forward, inverse))
secax.set_xlabel('CO$_2$ level (ppmv)')
secax.set_xticks(np.arange(400, 1400, 400))
fig.text(0.87, 0.73, '(h)', ha='center')
fig.tight_layout()
plt.show()
    
# savepath = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/Figures/Paper/Figure 1/'
# filename = savepath+'all_timeplots_V3.svg'
# fig.savefig(filename, format = 'svg', dpi=300, bbox_inches='tight') 