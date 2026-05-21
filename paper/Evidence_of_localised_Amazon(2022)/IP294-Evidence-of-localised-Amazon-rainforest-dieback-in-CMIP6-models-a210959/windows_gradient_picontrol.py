# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 12:01:16 2021

@author: Paul

Generate control run data relevant for determining abrupt shifts using the 
detection algorithm third criterion.
"""

import netCDF4 as nc
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.io import savemat

# Name of control run experiment
experiment = 'piControl'

################# Specify model here ################################## 
model = 'GFDL-ESM4'

################ Specify variant_id of model used ###################
variant_id = 'r1i1p1f1'

################# Specify variable of interest ######################
# var2 is the same as var but all lower case
var = 'tas'
var2 = 'tas'

################# Specify date range in file name #####################
date_range = '0000-1000'

# Is variable measured in per second?
PERSEC = False

####### Specify name of directory of data which you want to run algorithm on ############
region = 'Amazon'#'World'#

############ Specify name of directory of that stores analysis data #################
region2 = 'Amazon'#'World'#

################## Specify path to processed data directory #####################
path = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment+'/Processed_data/'+region+'/'

################## Specify path to analysis data directory #####################
path2 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment+'/Analysis_data/'+region2+'/'

# Determine window length for period of abrupt shift
wl = 15

# File name of the data
fname = path+var+'_'+model+'_'+experiment+'_'+variant_id+'_'+date_range+'_'+region+'.nc'

# Extract data
f = nc.Dataset(fname,'r')

# Extract longitudes and latitudes
longitude = f.variables['lon'][:]
latitude = f.variables['lat'][:]

# Extract variable data and convert to yearly units if applicable
if PERSEC:
    x = f.variables[var2][:,:,:]*3600*24*360
else:
    x = f.variables[var2][:,:,:]

# Close dataset
f.close()

# Determine dimensions of 3d array though use only first 150 years of data
nt = int(len(x[:,0,0]))
ny = int(len(x[0,:,0]))
nx = int(len(x[0,0,:]))

############ Specify path to land mask ######################
fname2 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/Masks/'+region+'.nc'

# Load in land mask and remove any ocean values
f2 = nc.Dataset(fname2,'r')
sftlf = f2.variables['sftlf'][:]
sftlf_mask = np.broadcast_to(sftlf, (nt, ny,nx))
x[np.where(sftlf_mask <= 0)] = np.nan  
x[np.where(~np.isfinite(sftlf_mask))] = np.nan
f2.close()

# Time parameters
t2 = np.linspace(0, wl-1, wl)

# Initialise array
control_grads = np.zeros((nt,ny,nx))

# Initialise counter
count = 0

# Create meshgrid of longitude and latitude
Lon, Lat = np.meshgrid(longitude, latitude)

# Loops over each gridpoint    
for i in range(nx):
    for j in range(ny):
        
        # If ocean all arrays are assigned a NaN 
        if not np.isfinite(x[0,j,i]) or all(x[:,j,i] == 0.0):
            control_grads[:,j,i] = np.nan
        
        else:

            # Loop over a sliding window
            for k in range(nt-wl):
                
                # Sliding window of variable
                x2 = x[k:wl+k,j,i]
                
                # Linear regression over sliding window
                regressor2 = LinearRegression()  
                regressor2.fit(t2.reshape(-1,1), x2)
                
                # From linear regression determine absolute change over sliding window
                control_grads[k,j,i] = regressor2.coef_

# Create dictionary of data to save                
mdic = {'Lon': Lon, 'Lat': Lat, 'control_grads': control_grads}

############# Save data as a matifle - Specify filename #############
savemat(path2+var+'_'+model+'_wl'+str(wl)+'_data_'+experiment+'_'+region2+'.mat', mdic)
