# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:45:14 2021

@author: impy2

Extracting the seasonal cycle anomalies from the data to remove the seasonal
cyle from considerations
"""
from pathlib import Path
import pandas as pd
import netCDF4 as nc
import numpy as np
from auto_corr import auto_corr
from scipy import stats
from scipy.io import loadmat
from scipy.io import savemat

################# Specify the experiment and models ################

experiment = '1pctCO2'
variant_id = 'r1i1p1f1'

################# Specify variable of interest ######################
var = 'tas' # or treeFrac
var2 = 'tas' # or treefrac
tip_var = 'cVeg'

############## Specify units of variable ###########################
units = '$kg/m^2$'

################# Specify model here ################################## 
model = 'TaiESM1'

############# Specify the longitude and latitude of the tipping point you wish to extract #########
latitude = -0
longitude = -60

######################## choose how much to calculate ########################
monthly_autocorr = True
diff_autocorr = True

if model == 'EC-Earth3-Veg':

    date_range = '1850-2000'
    total_years = 151
    
elif model == 'GFDL-ESM4':
    #modified to fit date range on downloaded data 
    date_range = '1850-1999'
    total_years = 150
    
elif model == 'MPI-ESM1-2-LR':
    
    date_range = '1850-2014'
    total_years = 165

elif model == 'NorCPM1':
    
    date_range = '1850-2013'
    total_years = 164

elif model == 'TaiESM1':
    
    date_range = '1850-2000'
    total_years = 150

elif model == 'SAM0-UNICON':
    
    date_range = '1850-1999'
    total_years = 150
    
elif model == 'UKESM1-0-LL':
    date_range = '1850-2000'
    total_years = 150
    
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
path2 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+tip_var+'/'+experiment+'/Analysis_data/'+region2+'/'
path3 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment+'/Statistical_data_monthly/'+region+'/'

# File name of interpolated data
fname = path+var+'_'+model+'_'+experiment+'_'+variant_id+'_'+date_range+'_mon_'+region+'.nc'
fname = Path(fname)

# # Load in data
f = nc.Dataset(fname,'r')
mat = loadmat(path2+tip_var+'_'+model+'_as_grads_wl'+str(wl)+'_data_'+experiment+'_'+region2+'.mat')

#extract indicies of lon and lat in the matrices 
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

data = x[:, lat_point, lon_point]
data = data.flatten()
t = np.arange(1, data.size+1)

#statistical analysis of resulting data 
co2_start = 284.3186666723341
co2 = np.zeros(data.size)
for k in range(data.size):
    co2[k] = co2_start*(1.01**k)
    
#extracting the CO2 conc and time point at which abrupt shift is detected 
co2_tip = mat['co2_tip']
co2_at_tip = co2_tip[lat_point, lon_point]
tip_index = np.where(co2==co2_at_tip)
time_at_tip = t[tip_index[1]] 

dataf = data

if monthly_autocorr == True:
    ####################### perform autocorrelation on the monthly data ##################
    window = 360
    ttot = np.size(t)
    mwind = ttot-window
    ar1_1 = np.zeros(mwind)
    stdev=np.zeros(mwind)
    ar1=np.zeros(mwind)
    skew=np.zeros(mwind)
    t1=t[0:mwind]
    stdev_d=np.zeros(mwind)
    ar1_d=np.zeros(mwind)
    skew_d=np.zeros(mwind)
    var1 = np.zeros(mwind)
    var_d = np.zeros(mwind)
    
    for m in range (0, mwind):
        m1 = m
        m2 = m+window
        veg = dataf[m1:m2]
        time = t[m1:m2]
        ####################### remove seasonal cycle ##############################
        #create arrays of data point indicies for each month
        month_points = np.zeros((12, 30))
        monthly_data = np.zeros((12,30))
        mon_avr = np.zeros(12)
        data2 = np.zeros(360)
        
        # extract monthly data points from the data
        for j in range (0,12):
            month_points[j, :] = np.arange(j, veg.size, 12)
            for i in range (0,30):
                monthly_data[j,i] = veg[int(month_points[j,i])]
                
            # find the average of each month
            mon_avr[j] = np.average(monthly_data[j,:])
        
         #subtract the average from each data point
        for j in range (0, 12):
            for i in range (0, 30):
                       data2[int(month_points[j,i])] = veg[int(month_points[j,i])] - mon_avr[j]
        
        b0, a0, r0, p0, e0 = stats.linregress(time,data2)
        dx=data2-(a0+b0*time) 
        stdev[m]=np.std(dx)
        var1[m] = stdev[m]**2
        ar1_x, dar1_x = auto_corr(dx,1)
        ar1_1[m]=ar1_x 
        skew_x=stats.skew(dx, bias=True)
        skew_d[m]=skew_x
        
if diff_autocorr == True:    
    #################### calculate the change in vegetation carbon ###############
    data = pd.Series(data)
    veg_flux = data.diff()
    
    #################### cut off data at abrupt shift ############################
    veg_flux_f = veg_flux
    flux_t = np.arange(1, veg_flux_f.size+1)
     
    #################### calculate variance and autocorrelation ##################
    window = 360
    ttot = np.size(flux_t)
    mwind = ttot-window
    ar1 = np.zeros(mwind)
    stdev=np.zeros(mwind)
    ar1_2=np.zeros(mwind)
    skew=np.zeros(mwind)
    t1=flux_t[0:mwind]
    var2 = np.zeros(mwind)
    
    for m in range (0, mwind):
        m1 = m
        m2 = m+window
        veg = veg_flux_f[m1:m2]
        veg = veg.to_numpy()
        time = flux_t[m1:m2]
        ####################### remove seasonal cycle ##############################
        #create arrays of data point indicies for each month
        month_points = np.zeros((12, 30))
        monthly_data = np.zeros((12,30))
        mon_avr = np.zeros(12)
        data3 = np.zeros(360)
        
        # extract monthly data points from the data
        for j in range (0,12):
            month_points[j, :] = np.arange(j, veg.size, 12)
            for i in range (0,30):
                monthly_data[j,i] = veg[int(month_points[j,i])]
                
            # find the average of each month
            mon_avr[j] = np.average(monthly_data[j,:])
        
         #subtract the average from each data point
        for j in range (0, 12):
            for i in range (0, 30):
                       data3[int(month_points[j,i])] = veg[int(month_points[j,i])] - mon_avr[j]
        
        #################### signal detrend ##################################
        b0, a0, r0, p0, e0 = stats.linregress(time,data3)
        dx=data3-(a0+b0*time) 
        stdev[m]=np.std(dx)
        var2[m] = stdev[m]**2
        ar1_x, dar1_x = auto_corr(dx,1)
        ar1_2[m]=ar1_x 
        skew_x=stats.skew(dx, bias=True)
        skew[m]=skew_x
 
################ save as mat files #######################################    
norm_var1 = var1/var1[1]
norm_var2 = var2/var2[1]

# Create dictionary of data to save                
mdic = {'time_tas':dataf, 'ar1': ar1_1, 'flux_ar1': ar1_2, 'var': norm_var1, 'flux_var': norm_var2, 'skew':skew}

############# Save data as a matifle - Specify filename #############
# savemat(path3+var+'_'+model+'_'+experiment+'_'+region2+'_lat'+str(latitude)+'_lon'+str(longitude)+'_statistical_data.mat', mdic)