# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 15:00:25 2021

@author: impy2

Plot collating all statistical analysis for generic EWS analysis
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy.io import loadmat
plt.rcParams.update({'font.size': 10})

################# Specify the experiment and models ################

experiment = '1pctCO2'
variant_id = 'r1i1p1f1'

################# Specify variable of interest ######################
# var2 is the same as var but all lower case
v_var = 'cVeg' # or treeFrac
t_var = 'tas' # or treefrac
var2 = 'cveg'

############## Specify units of variable ###########################
units = '$kg/m^2$'

################# Specify models here ################################## 
models = ['GFDL-ESM4', 'NorCPM1', 'TaiESM1'] 
##panel labels for subplots
alphabet = ['d', 'e', 'f']

############# Specify the longitude and latitude of the tipping point you wish to extract #########
latitudes = [-5, 0, 0]
longitudes = [-65, -60, -60]

###################### window size ###########################################
window = 360

###################### specify region #######################################
region = 'Amazon'

###################### specify variable ######################################
var1 = 'tas'
var2 = 'cveg'

#################### load in mat files #######################################
path1 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var1+'/'+experiment+'/Statistical_data_monthly/'+region+'/'
path2 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var2+'/'+experiment+'/Statistical_data_monthly/'+region+'/'

## initialise figure
fig, axes = plt.subplots(3, 1, figsize=(4.1, 4.9), sharex=True)#, sharey=True)

##loop through models
for i in range (0,3):
    model = models[i]
    latitude = latitudes[i]
    longitude = longitudes[i]
    tmat = loadmat(path1+var1+'_'+model+'_'+experiment+'_'+region+'_lat'+str(latitude)+'_lon'+str(longitude)+'_statistical_data.mat')
    cmat = loadmat(path2+var2+'_'+model+'_'+experiment+'_'+region+'_lat'+str(latitude)+'_lon'+str(longitude)+'_statistical_data.mat')
    rmat = loadmat(path2+var2+'_'+model+'_'+experiment+'_'+region+'_lat'+str(latitude)+'_lon'+str(longitude)+'_linreg.mat')
    ######################### define variables ##################################
    cveg_data = cmat['cveg']
    cflux_data = cmat['cflux']
    cveg_ar1 = cmat['ar1']
    tmat['time_tas']

    t = np.arange(1, np.size(cmat['cveg'])+1)
    t_adj = np.arange(1, np.size(cmat['ar1'])+1)
    t_adj = (t_adj+window)/12
  
    ####### variance variables #######
    y1=cmat['var']/cmat['var'][0,1]
    
    ##### autocorrelation variables #####
    y2 = cmat['ar1'].flatten()
    
    # Create array of CO2 for 1pctCO2 run
    co2_start = 284.3186666723341
    co2 = np.zeros(int(y1.size/12))
    for k in range(int(y1.size/12)):
        co2[k] = co2_start*(1.01**k)
        
    time = np.arange(1, np.size(y1)+1)
    
    ## create the secondary axis data functions ####
    # time to co2
    def forward(x):
        return co2_start*(1.01**x)
    
    # co2 to time
    def inverse(x):
        return np.log(x/co2_start)/np.log(1.01)
   
    ####### plot figure ########
    plt.xlabel('Time (model years)')
    axes[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axes[i].text(0.92, 0.82, '('+alphabet[i]+')', transform=axes[i].transAxes)
    axes[i].plot(t_adj[0:-12], y1.flatten()[0:-12], label= 'Variance', color='green', linestyle='dotted')
    ax2 = axes[i].twinx()
    ax2.plot(t_adj, y2.flatten()[0:t_adj.size], label = 'AR1', color='green', linestyle='dashed')
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    if i ==0:
        axes[i].set_ylim(0,20)
        lines, labels = axes[i].get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc='upper left') 
    
plt.xticks(np.arange(0, 160,20))
plt.xlim(0,150)
secax = axes[0].secondary_xaxis('top', functions=(forward, inverse))
plt.draw()
secax.set_xticklabels(secax.get_xticks().astype(int), rotation =45)
secax.set_xlabel('CO$_2$ level (ppmv)')

fig.text(-0.01, 0.5, 'Variance (normalised)', va='center', rotation='vertical')
fig.text(1.03, 0.5, 'AR1', va='center', rotation='vertical')

# plt.legend(ncol=2,loc='best', bbox_to_anchor=(1, -0.5, -0., 0.))
# savepath = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/Figures/Paper/Supplementary/EWS/'
# filename = savepath+'EWS.svg'
# fig.savefig(filename, format = 'svg', dpi=300, bbox_inches='tight') 
    
 