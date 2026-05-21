# -*- coding: utf-8 -*-
"""
Created on Sun May 16 09:07:13 2021

Create a plot of the yearly maximum mean month temp 

@author: impy2
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})
from scipy.io import loadmat
from matplotlib.ticker import FormatStrFormatter

################# Specify the experiment and models ################

experiment = '1pctCO2'
variant_id = 'r1i1p1f1'

################# Specify variable of interest ######################
var = 'tas'

############## Specify units of variable ###########################
units = '$kg/m^2$'

################# Specify model here ################################## 
model1 = 'GFDL-ESM4'
model2 = 'NorCPM1'
model3 = 'TaiESM1'

alphabet = ['d', 'e', 'f']

############# Specify the longitude and latitude of the tipping point you wish to extract #########
latitude1 = -5
longitude1 = -65

latitude2 = 0
longitude2 = -60

latitude3 = 0
longitude3 = -60

###################### window size ###########################################
window = 360

###################### specify region #######################################
region = 'Amazon'

###################### specify variable ######################################
var1 = 'tas'

#################### load in mat files #######################################
path1 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var1+'/'+experiment+'/Statistical_data_monthly/'+region+'/'

tmat1 = loadmat(path1+var1+'_'+model1+'_'+experiment+'_'+region+'_lat'+str(latitude1)+'_lon'+str(longitude1)+'_statistical_data.mat')
tmat2 = loadmat(path1+var1+'_'+model2+'_'+experiment+'_'+region+'_lat'+str(latitude2)+'_lon'+str(longitude2)+'_statistical_data.mat')
tmat3 = loadmat(path1+var1+'_'+model3+'_'+experiment+'_'+region+'_lat'+str(latitude3)+'_lon'+str(longitude3)+'_statistical_data.mat')

T1 = tmat1['time_tas']
T1 = T1.flatten()

T2 = tmat2['time_tas']
T2 = T2.flatten()

T3 = tmat3['time_tas']
T3 = T3.flatten()

#### separate data into years #####
sets = [T1, T2, T3]
lats = [latitude1, latitude2, latitude3]
lons = [longitude1, longitude2, longitude3]
models = [model1, model2, model3]

##specify time at which abrupt shift occured for each
time_at_tip = [91, 107, 90]

## initialise figure
fig, axes = plt.subplots(3, 1, figsize=(4.1, 4.9), sharex=True)

#loop over each model to create subplots
for j in range (0, 3):
   
    s = int(np.size(sets[j])/12)
    year_points = np.zeros((s, 12))
    max_y = np.zeros(s)
    mean_y = np.zeros(s)
    diff_y = np.zeros(s)
    min_y = np.zeros(s)
    amp = np.zeros(s)
    run_amp = np.zeros(s-10)

    for i in range (0, s):
        year_points[i, :] = np.arange(i*12, (i*12)+12)
    
    #extract data per year and then find the max temp in each yr 
    for i in range (0,s):
        points = year_points[i,:]
        d_y = sets[j][int(points[0]):int(points[-1])+1]
        max_y[i] = np.max(d_y)
        mean_y[i] = np.average(d_y)
        diff_y[i] = max_y[i] - mean_y[i]
        #add in calc for the full seasonal cycle amplitude
        min_y[i] = np.min(d_y)
        amp[i] = max_y[i] - min_y[i]
        
        ### calculate the running average ###
    for k in range (0, (np.size(amp)-10)):
        x1 = amp[k:k+10]
        run_amp[k] = np.average(x1)
        
    # Create array of CO2 for 1pctCO2 run
    co2_start = 284.3186666723341
    co2 = np.zeros(int(amp.size/12))
    for k in range(int(amp.size/12)):
        co2[k] = co2_start*(1.01**k)
        
    time = np.arange(1, np.size(amp)+1)
    time2 = (np.arange(1, np.size(run_amp)+1))+10
    
    ## create the secondary axis data functions ####
    # time to co2
    def forward(x):
        return co2_start*(1.01**x)
    
    # co2 to time
    def inverse(x):
        return np.log(x/co2_start)/np.log(1.01)


    if lats[j]>0:
        lat_str = str(int(np.linalg.norm(lats[j])))+u'\N{DEGREE SIGN}N'
    else:
        lat_str = str(int(np.linalg.norm(lats[j])))+u'\N{DEGREE SIGN}S'
        
    lon_str = str(int(np.linalg.norm(lons[j])))+u'\N{DEGREE SIGN}W'

    if j==1:
        axes[j].plot(time2[0:138], run_amp[0:138], color='goldenrod')
        axes[j].plot(time[0:148], amp[0:148], color='goldenrod', linestyle='dotted') 
    else:
        axes[j].plot(time2, run_amp, color='goldenrod')
        axes[j].plot(time, amp, color='goldenrod', linestyle='dotted')
        
    axes[j].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axes[j].vlines(time_at_tip[j], min(amp), max(amp), colors='red', linestyles='dashed')
    
    axes[j].text(0.92, 0.08, '('+alphabet[j]+')', transform=axes[j].transAxes)
    axes[j].set_yticklabels(axes[j].get_yticks().astype(int))
    
    if j == 0:
        secax = axes[0].secondary_xaxis('top', functions=(forward, inverse))
        plt.draw()
        secax.set_xticklabels(secax.get_xticks().astype(int), rotation =45)
        secax.set_xlabel('CO$_2$ level (ppmv)')
    if models[j] == 'GFDL-ESM4':
        axes[j].set_yticks(np.arange(0, 16, 4 ))
    if models[j] == 'NorCPM1':
        axes[j].set_yticks(np.arange(1, 7, 2))
        
plt.xlabel('Time (model years)')
fig.text(-0.01, 0.5, 'T Seasonal Cycle Amplitude (K)', va='center', rotation='vertical')
plt.xticks(np.arange(0,160,20))
plt.show()
    
# savepath = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/Figures/Paper/Figure 3/'
# filename = savepath+'temp_plot(3).svg'
# fig.savefig(filename, format = 'svg', dpi=300, bbox_inches='tight') 






