# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 12:55:01 2021

@author: impy2

Script to create time-series plots for specifc tipping points  as in figure 3
"""

from pathlib import Path
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from matplotlib.ticker import FormatStrFormatter
plt.rcParams.update({'font.size': 10})

################# Specify the experiment and models ################

experiment = '1pctCO2'
variant_id = 'r1i1p1f1'

################# Specify variable of interest ######################
# var2 is the same as var but all lower case
var = 'cVeg' # or treeFrac
var2 = 'cveg' # or treefrac

############## Specify units of variable ###########################
units = '$kgC m^{-2}$'

## define figure 
fig, axes = plt.subplots(3, 1, figsize=(4.1, 4.9), sharex=True)
models = ['GFDL-ESM4', 'NorCPM1', 'TaiESM1']
alphabet = ['a', 'b', 'c']

lats = [-5, 0, 0]
lons = [-65, -60, -60]

#loop for each model
for i in range (0,3):
    model = models[i]
    latitude = lats[i]
    longitude = lons[i]
    
    if model == 'EC-Earth3-Veg':
    
        date_range = '1850-2000'
        
    elif model == 'GFDL-ESM4':
        #modified to fit date range on downloaded data 
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
    
    # Create array of CO2 for 1pctCO2 run
    co2_start = 284.3186666723341
    co2 = np.zeros(int(data.size/12))
    for k in range(int(data.size/12)):
        co2[k] = co2_start*(1.01**k)
    
    # #extracting the CO2 conc and time point at which abrupt shift is detected 
    co2_tip = mat['co2_tip']
    co2_at_tip = co2_tip[lat_point, lon_point]
    tip_index = np.where(co2==co2_at_tip)
    time_at_tip = t[tip_index[1]]     
     
    data_yr = np.zeros(int(np.size(data)/12))
    ## average the monthly data for each year
    for k in range (0, (int((np.size(data)/12)))):
        x1 = data[(k*12):(k*12)+12]
        data_yr[k] = np.average(x1)                
    time = np.arange(0, data_yr.size,1)
                         
    ############## plot the time series data ############
    if latitude>0:
        lat_str = str(int(np.linalg.norm(latitude)))+u'\N{DEGREE SIGN}N'
    else:
        lat_str = str(int(np.linalg.norm(latitude)))+u'\N{DEGREE SIGN}S'
        
    lon_str = str(int(np.linalg.norm(longitude)))+u'\N{DEGREE SIGN}W'
    
    ## create the secondary axis data functions ####
    # time to co2
    def forward(x):
        return co2_start*(1.01**x)
    
    # co2 to time
    def inverse(x):
        return np.log(x/co2_start)/np.log(1.01)
    
    if model == 'NorCPM1':
        axes[i].plot(time[0:148], data_yr[0:148], 'g')
    else:
        axes[i].plot(time, data_yr, 'g')
        
    axes[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axes[i].vlines(time_at_tip, min(data), max(data), colors='red', linestyles='dashed')
    axes[i].text(0.92, 0.82, '('+alphabet[i]+')', transform=axes[i].transAxes)
     
    #draw secondary x axis
    if i == 0:
        secax = axes[0].secondary_xaxis('top', functions=(forward, inverse))
        plt.draw()
        secax.set_xticklabels(secax.get_xticks().astype(int), rotation =45)
        secax.set_xlabel('CO$_2$ level (ppmv)')   
    
    axes[i].set_yticklabels(axes[i].get_yticks().astype(int))
    axes[i].set_title(model, loc='left', x=0.03, y=0.05)
    
plt.xlabel('Time (model years)')
plt.xticks(np.arange(0,160,20))
fig.text(0, 0.5, 'Vegetation Carbon' +' ('+units+')', va='center', rotation='vertical')
plt.show()

# savepath = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/Figures/Paper/Figure 3/'
# filename = savepath+'monthly_cveg_plots(4).svg'
# fig.savefig(filename, format = 'svg', dpi=300, bbox_inches='tight') 


    
