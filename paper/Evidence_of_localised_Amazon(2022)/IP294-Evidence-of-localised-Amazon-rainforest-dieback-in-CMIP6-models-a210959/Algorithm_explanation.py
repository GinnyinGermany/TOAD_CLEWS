# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 11:53:51 2022

@author: impy2
"""
from pathlib import Path
# import iris
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from scipy.io import loadmat
plt.rcParams.update({'font.size': 10})

################# Specify the experiment and models ################
experiment = '1pctCO2'
experiment2 = 'piControl'
variant_id = 'r1i1p1f1'

################# Specify variable of interest ######################
var = 'cVeg' # or treeFrac
var2 = 'cveg' # or treefrac

############## Specify units of variable ###########################
units = '$kgC m^{-2}$'
units2 = '$kgC m^{-2} yr^{-1}$'

################# Specify model here ################################## 
model = 'EC-Earth3-Veg'

############# Specify the longitude and latitude of the tipping point you wish to extract #########
latitude = -10
longitude = -68

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

#set the threshold number of standard deviations to use for criteria
std_threshold = 3

################## Specify path to processed data directory #####################
path = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment+'/Processed_data/'+region+'/'
path2 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment+'/Analysis_data/'+region2+'/'

# File name of interpolated data
fname = path+var+'_'+model+'_'+experiment+'_'+variant_id+'_'+date_range+'_'+region+'.nc'
fname = Path(fname)

## Load in data
f = nc.Dataset(fname,'r')

## load in mat file created by abrupt_shift_detection_statistics 
mat = loadmat(path2+var+'_'+model+'_as_grads_wl'+str(wl)+'_data_'+experiment+'_'+region2+'.mat')

##extract indicies of lon and lat in the matrices 
lon = np.array(f.variables['lon'])
lon_point = np.where(lon==longitude)

lat = np.array(f.variables['lat'])
lat_point = np.where(lat==latitude)

## Extract data and if measured in per second convert to per year
if PERSEC:
    x = f.variables[var2][:]*3600*24*360
else:
    x = f.variables[var2][:]

# Close dataset
f.close()

### extract the data as a time series ####

data = x[:, lat_point, lon_point]
data = data.flatten()
## remove last data point
data=data[0:-1]
t = np.arange(1, data.size+1)

## Create array of CO2 for 1pctCO2 run
co2_start = 284.3186666723341
co2 = np.zeros(data.size)
for k in range(data.size):
    co2[k] = co2_start*(1.01**k)

## extract the CO2 conc and time point at which abrupt shift is detected 
co2_tip = mat['co2_tip']
co2_at_tip = co2_tip[lat_point, lon_point]
tip_index = np.where(co2==co2_at_tip)
time_at_tip = t[tip_index[1]]                     

## load in mat file created by abrupt_shift_detection_statistics                       
path4 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment2+'/Analysis_data/'+region2+'/'
mat2 = loadmat(path4+var+'_'+model+'_wl'+str(wl)+'_data_'+experiment2+'_'+region2+'.mat')
    
control_grads = mat2['control_grads']
control_grad_std = np.nanstd(control_grads[-400:-wl,:,:], axis=0)
std_1 = control_grad_std[lon_point,lat_point]

#### take gradient of the line of best fit over the 15 year window ####
t2 = np.linspace(0, wl-1, wl)
as_grads = np.zeros(np.size(t)-wl)
#loop over sliding window
for k in range(np.size(t)-wl):
     # Sliding window of variable
    data2 = data[k:wl+k]
    # Linear regression over sliding window
    regressor = LinearRegression()
    regressor.fit(t2.reshape(-1,1), data2)
    # From linear regression determine absolute change over sliding window
    as_grads[k] = regressor.coef_

dydx = np.gradient(data, t)

############## plot the time series data ############
if latitude>0:
    lat_str = str(int(np.linalg.norm(latitude)))+u'\N{DEGREE SIGN}N'
else:
    lat_str = str(int(np.linalg.norm(latitude)))+u'\N{DEGREE SIGN}S'
    
lon_str = str(int(np.linalg.norm(longitude)))+u'\N{DEGREE SIGN}W'

########### create figure ############
fig, axes = plt.subplots(2,1, figsize=(4.1, 3.5), sharex=True)

fig.suptitle(model)
plt.tight_layout()
axes[0].plot(t, data, 'g')
axes[0].set_ylabel(var+' ('+units+')')
axes[0].set_ylim([0,30])

axes[1].hlines(std_1*std_threshold, np.min(t), np.max(t), linestyle='dashed')
axes[1].hlines(-std_1*std_threshold, np.min(t), np.max(t), linestyle='dashed')
axes[1].plot(t[7:-8], as_grads, 'g')
axes[1].set_ylabel('d('+var+')/dt \n' +units2)
axes[1].set_ylim([-1.1,0.5])
axes[1].set_yticks([-1, -0.5, 0, 0.5])

axes[0].text(0.92, 0.82, '(a)', transform=axes[0].transAxes)
axes[1].text(0.92, 0.82, '(b)', transform=axes[1].transAxes)


plt.xlabel('Time (model years)')
plt.xticks(np.arange(0, np.size(t)+10, step=20)) 
fig.tight_layout()

# savepath = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/Figures/Paper/Supplementary/'
# filename = savepath+'Algorithm_std_'+model+'_'+str(latitude)+'_'+str(longitude)+'.svg'
# fig.savefig(filename, format = 'svg', dpi=300, bbox_inches='tight') 
