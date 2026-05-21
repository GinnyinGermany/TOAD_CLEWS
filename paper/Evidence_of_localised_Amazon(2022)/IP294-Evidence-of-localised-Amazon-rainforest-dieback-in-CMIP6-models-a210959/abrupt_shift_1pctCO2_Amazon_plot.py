# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 11:13:29 2021

@author: Paul

Plots for simplified algorithm for detecting abrupt shift: abrupt shift must
contribute to at least X % of overall shift and have at least a Y % change over
15 years of time series
"""

import numpy as np
import matplotlib.pyplot as plt

# plt.rcParams.update({'font.size': 23})

import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker
from scipy.io import loadmat
import matplotlib.patches as mpatches
from scipy.io import savemat


################# Specify variable of interest ######################
var = 'cVeg'

################# Specify model here ################################## 
model = 'TaiESM1'

######## Determine if control statistics are to be used for plotting #########
######## detected abrupt shifts ##############################################
CONTROL = True

# Specify experiments
experiment = '1pctCO2'
experiment2 = 'piControl'

####### Specify name of directory of data which you want to run algorithm on ############
region = 'Amazon'#'World'#

############ Specify name of directory of that stores analysis data #################
region2 = 'Amazon'#'World'#

# Longitude, latitude min, max, step and figure size for region of interest
if region == 'World':
    lonmin = -180-0.5
    lonmax = 180+0.5
    latmin = -60-0.5
    latmax = 80+0.5
    lonstep = 10
    latstep = 10
    figsizes = (19, 9.5)
elif region == 'Amazon':
    lonmin = -83-0.5
    lonmax = -35+0.5
    latmin = -20-0.5
    latmax = 13+0.5
    lonstep = 10
    latstep = 10
    figsizes = (6, 5.75)

# Window length used to calculate abrupt shift over
wl = 15

################## Specify path to analysis data directory #####################
path = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment+'/Analysis_data/'+region2+'/'
path2 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment2+'/Analysis_data/'+region2+'/'
path3 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment+'/tip_points_data/'+region2+'/'

################## Specify filename of saved analysis data and load in #####################
mat = loadmat(path+var+'_'+model+'_as_grads_wl'+str(wl)+'_data_'+experiment+'_'+region2+'.mat')
mat2 = loadmat(path2+var+'_'+model+'_wl'+str(wl)+'_data_'+experiment2+'_'+region2+'.mat')

# Extract data from dictionary
Lon, Lat, as_grad, as_grads, co2_tip, as_change, ovr_change = mat['Lon'], mat['Lat'], mat['as_grad'], mat['as_grads'], mat['co2_tip'], mat['as_change'], mat['ovr_change']

# Load control run data in if applicable
if CONTROL:
    ################## Specify path to analysis data directory #####################
    path2 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment2+'/Analysis_data/'+region2+'/'

    ################## Specify filename of saved analysis data and load in #####################
    mat2 = loadmat(path2+var+'_'+model+'_wl'+str(wl)+'_data_'+experiment2+'_'+region2+'.mat')

    control_grads = mat2['control_grads']
    control_grad_std = np.nanstd(control_grads[-400:-wl,:,:], axis=0)

    ############ Set minimum number of standard deviations from ############
    ############ zero for it to be classed as an abrupt shift ############
    std_threshold = 3

############ Set minimum threshold for change in sliding window to #############
############ be classed as an abrupt shift #####################################
abs_change_threshold = 2

############ Set minimum fractional contribution of abrupt shift to ##########
############ overall change ##################################################
as_frac_threshold = 0.25

# Marker and label size for colour scatter plots
if region2 == 'World':
    marker_size = 4
    proj = 'cyl'
    loc = 'lower left'
    marker_label_size = 6
else:
    marker_size = 40.0
    proj = 'merc'
    loc = 'upper right'
    marker_label_size = 2

# Labels of abrupt shift types and their colours
type_labels = ['T>0, AS>0', 'T<0, AS<0', 'T<0, AS>0', 'T>0, AS<0']
colours = ['g','r','b','tab:orange']

step = 0.5

# Colormap for abrupt shift change
v1  = [-10, -5,-3, -2, 2,3, 5, 10]#[-1.7, -1.4, -1.1,-0.8, 0.8, 1.1, 1.4, 1.7]#[-1000, -800, -600,-500, 500, 600, 800, 1000]#[-30, -20, -10, 10, 20, 30]#
cmap1_name = 'RdYlGn'
cmap1 = plt.cm.get_cmap(cmap1_name,len(v1)+1)
colors1 = list(cmap1(np.arange(len(v1)+1)))
cmap1 = colors.ListedColormap(colors1[1:-1], "")
cmap1.set_under(colors1[0])     # set over-color to last color of list 
cmap1.set_over(colors1[-1])     # set under-color to firstst color of list 
norm1 = colors.BoundaryNorm(v1, cmap1.N)

# Colormap for no. of std above control
v2 = [0, 1, 3, 5, 10, 50, 100]
cmap2_name = 'plasma'
cmap2 = plt.cm.get_cmap(cmap2_name,len(v2))
colors2 = list(cmap2(np.arange(len(v2))))
cmap2 = colors.ListedColormap(colors2[:-1], "")
cmap2.set_over(colors2[-1])     # set over-color to last color of list 
norm2 = colors.BoundaryNorm(v2, cmap2.N)

# Colormap for abrupt shift change
v3  = [-10, -5,-3, -2, 2,3, 5, 10]#[-30, -20, -10, 10, 20, 30]#
cmap3_name = 'RdYlGn'
cmap3 = plt.cm.get_cmap(cmap3_name,len(v3)+1)
colors3 = list(cmap3(np.arange(len(v3)+1)))
cmap3 = colors.ListedColormap(colors3[1:-1], "")
cmap3.set_under(colors3[0])     # set over-color to last color of list 
cmap3.set_over(colors3[-1])     # set under-color to firstst color of list 
norm3 = colors.BoundaryNorm(v3, cmap3.N)

# Colormap for CO2 tipping time
co2_start = 284.3186666723341
v4 = [co2_start]
v4.extend(np.arange(1.25, 4.75, step)*co2_start)
cmap4 = colors.ListedColormap(["darkred", "red", "darkorange", "gold", "darkgreen", "lime", "aqua"])
norm4 = colors.BoundaryNorm(v4, cmap4.N)

# Colormap for abrupt shift fraction
v5 = [0, 0.25, 0.5, 1, 2]
cmap5_name = 'plasma'
cmap5 = plt.cm.get_cmap(cmap5_name,len(v5))
colors5 = list(cmap5(np.arange(len(v5))))
cmap5 = colors.ListedColormap(colors5[:-1], "")
cmap5.set_over(colors5[-1])     # set over-color to last color of list 
norm5 = colors.BoundaryNorm(v5, cmap5.N)


if CONTROL:
    indxs = np.where((np.abs(as_grad)>std_threshold*control_grad_std)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
    indxs2 = ~((np.abs(as_grad)>std_threshold*control_grad_std)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
    indxs3 = ((ovr_change>0)&(as_change>0)&(np.abs(as_grad)>std_threshold*control_grad_std)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
    indxs4 = ((ovr_change<0)&(as_change<0)&(np.abs(as_grad)>std_threshold*control_grad_std)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
    indxs5 = ((ovr_change<0)&(as_change>0)&(np.abs(as_grad)>std_threshold*control_grad_std)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
    indxs6 = ((ovr_change>0)&(as_change<0)&(np.abs(as_grad)>std_threshold*control_grad_std)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
    
    AS_STD = np.ma.masked_where(indxs2, (np.abs(as_grad)/control_grad_std))
else:
    indxs = np.where((np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
    indxs2 = ~((np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
    indxs3 = ((ovr_change>0)&(as_change>0)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
    indxs4 = ((ovr_change<0)&(as_change<0)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
    indxs5 = ((ovr_change<0)&(as_change>0)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
    indxs6 = ((ovr_change>0)&(as_change<0)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))

AS_CHANGE = np.ma.masked_where(indxs2, as_change)
OVR_CHANGE = np.ma.masked_where(indxs2, ovr_change)
CO2_TIP = np.ma.masked_where(indxs2, co2_tip)
AS_FRAC = np.ma.masked_where(indxs2, np.abs(as_change/(ovr_change)))
AS_FRAC = np.ma.masked_where(indxs2, np.abs(as_change/(ovr_change)))

###################### save information in a mat dataset #####################
mdic = {'tips': AS_CHANGE}
savemat(path3+var+'_'+model+'_'+experiment+'_'+region2+'_as_change_data.mat', mdic)

# Find norm error to coordinates points
a1 = abs(Lat[:,0]-(latmin+0.5))
a2 = abs(Lat[:,0]-(latmax-0.5))
a3 = abs(Lon[0,:]-(lonmin+0.5))
a4 = abs(Lon[0,:]-(lonmax-0.5))

# Find array indicies of nearest point to coordinates 
index_1 = a1.argmin()
index_2 = a2.argmin()
index_3 = a3.argmin()
index_4 = a4.argmin()

##################### Plotting figures - specify file names ##################

proj = ccrs.Mercator(central_longitude=-59, min_latitude=latmin, max_latitude=latmax)
# save_path = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/Figures/'+var+'/'+experiment+'/'+region2+'/Gradient_detection/'
save_path = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/Figures/Paper/Supplementary/Thresholds'

# # # ### plot abrubt shift change ###

fig, ax = plt.subplots(subplot_kw=dict(projection=proj),figsize=figsizes)
ax.set_extent([lonmin ,lonmax, latmin, latmax], crs=ccrs.PlateCarree())   
ax.add_feature(cfeature.COASTLINE, zorder=10)
ax.set_xticks(np.arange(-80,-30,lonstep),crs=ccrs.PlateCarree())
ax.set_xticklabels(np.arange(-80,-30,lonstep))
ax.set_yticks(np.arange(-20,20,latstep),crs=ccrs.PlateCarree())
ax.set_yticklabels(np.arange(-20,20,latstep))
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.grid(linewidth=0.5, color='black', linestyle='--')
im = ax.imshow(AS_CHANGE,cmap=cmap1,norm=norm1,transform=ccrs.PlateCarree(),extent=[lonmin ,lonmax, latmax, latmin])

cbar = fig.colorbar(im, pad=0.08, orientation='horizontal', spacing="proportional",extend='both')
cbar.set_label('Abrupt shift absolute change ($kg C/m^2$)')
plt.title(model)
fig.tight_layout()
# if CONTROL:
#     fig.savefig(save_path+'Control/'+var+'_as_change_asthresh'+str(abs_change_threshold)+'_fracthresh'+str(as_frac_threshold).replace('.', '')+'_stdthresh'+str(std_threshold)+'_'+experiment+'_'+region2+'_'+model+'.png', bbox_inches='tight')
# else:
#     fig.savefig(save_path+var+'_as_change_asthresh'+str(abs_change_threshold)+'_fracthresh'+str(as_frac_threshold).replace('.', '')+'_'+experiment+'_'+region2+'_'+model+'.png', bbox_inches='tight')

##### Plot overall change #####

fig, ax = plt.subplots(figsize=figsizes)
fig, ax = plt.subplots(subplot_kw=dict(projection=proj),figsize=figsizes)
ax.set_extent([lonmin ,lonmax, latmin, latmax], crs=ccrs.PlateCarree())   
ax.add_feature(cfeature.COASTLINE, zorder=10)
ax.set_xticks(np.arange(-80,-30,lonstep),crs=ccrs.PlateCarree())
ax.set_xticklabels(np.arange(-80,-30,lonstep))
ax.set_yticks(np.arange(-20,20,latstep),crs=ccrs.PlateCarree())
ax.set_yticklabels(np.arange(-20,20,latstep))
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.grid(linewidth=0.5, color='black', linestyle='--')
im = ax.imshow(OVR_CHANGE,cmap=cmap3,norm=norm3,transform=ccrs.PlateCarree(),extent=[lonmin ,lonmax, latmax, latmin])

cbar = fig.colorbar(im, pad=0.08, orientation='horizontal', spacing="proportional",extend='both')
cbar.set_label('Overall absolute change ($kg C/m^2$)')
plt.title(model)
fig.tight_layout()
# if CONTROL:
#     fig.savefig(save_path+'Control/'+var+'_ovr_change_asthresh'+str(abs_change_threshold)+'_fracthresh'+str(as_frac_threshold).replace('.', '')+'_stdthresh'+str(std_threshold)+'_'+experiment+'_'+region2+'_'+model+'.png', bbox_inches='tight')
# else:
#     fig.savefig(save_path+var+'_ovr_change_asthresh'+str(abs_change_threshold)+'_fracthresh'+str(as_frac_threshold).replace('.', '')+'_'+experiment+'_'+region2+'_'+model+'.png', bbox_inches='tight')

####### Plot CO2 level at abrupt shift ########

fig, ax = plt.subplots(figsize=figsizes)
fig, ax = plt.subplots(subplot_kw=dict(projection=proj),figsize=figsizes)
ax.set_extent([lonmin ,lonmax, latmin, latmax], crs=ccrs.PlateCarree())   
ax.add_feature(cfeature.COASTLINE, zorder=10)
ax.set_xticks(np.arange(-80,-30,lonstep),crs=ccrs.PlateCarree())
ax.set_xticklabels(np.arange(-80,-30,lonstep))
ax.set_yticks(np.arange(-20,20,latstep),crs=ccrs.PlateCarree())
ax.set_yticklabels(np.arange(-20,20,latstep))
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.grid(linewidth=0.5, color='black', linestyle='--')
im = ax.imshow(CO2_TIP,cmap=cmap4,norm=norm4,transform=ccrs.PlateCarree(),extent=[lonmin ,lonmax, latmax, latmin])

cbar = fig.colorbar(im, ticks=np.arange(1, 4.5, step)*co2_start, pad=0.08, orientation='horizontal', spacing="proportional")
cbar.set_label('Factor of $CO_2$ increase')
cbar.ax.set_xticklabels(['1', '1.5', '2', '2.5', '3', '3.5', '4',])
plt.title(model)
fig.tight_layout()

# if CONTROL:
#     fig.savefig(save_path+'Control/'+var+'_as_co2_asthresh'+str(abs_change_threshold)+'_fracthresh'+str(as_frac_threshold).replace('.', '')+'_stdthresh'+str(std_threshold)+'_'+experiment+'_'+region2+'_'+model+'.png', bbox_inches='tight')
# else:
#     fig.savefig(save_path+var+'_as_co2_asthresh'+str(abs_change_threshold)+'_fracthresh'+str(as_frac_threshold).replace('.', '')+'_'+experiment+'_'+region2+'_'+model+'.png', bbox_inches='tight')


####### Plot fraction abrupt shift contributes to overall change #######

fig, ax = plt.subplots(figsize=figsizes)
fig, ax = plt.subplots(subplot_kw=dict(projection=proj),figsize=figsizes)
ax.set_extent([lonmin ,lonmax, latmin, latmax], crs=ccrs.PlateCarree())   
ax.add_feature(cfeature.COASTLINE, zorder=10)
ax.set_xticks(np.arange(-80,-30,lonstep),crs=ccrs.PlateCarree())
ax.set_xticklabels(np.arange(-80,-30,lonstep))
ax.set_yticks(np.arange(-20,20,latstep),crs=ccrs.PlateCarree())
ax.set_yticklabels(np.arange(-20,20,latstep))
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.grid(linewidth=0.5, color='black', linestyle='--')
im = ax.imshow(AS_FRAC,cmap=cmap5,norm=norm5,transform=ccrs.PlateCarree(),extent=[lonmin ,lonmax, latmax, latmin])

cbar = fig.colorbar(im, pad=0.08, orientation='horizontal', spacing="proportional",extend='both')
cbar.set_label('Fraction of overall change')
plt.title(model)
fig.tight_layout()

# if CONTROL:
#     fig.savefig(save_path+'Control/'+var+'_as_frac_asthresh'+str(abs_change_threshold)+'_fracthresh'+str(as_frac_threshold).replace('.', '')+'_stdthresh'+str(std_threshold)+'_'+experiment+'_'+region2+'_'+model+'.png', bbox_inches='tight')
# else:
#     fig.savefig(save_path+var+'_as_frac_asthresh'+str(abs_change_threshold)+'_fracthresh'+str(as_frac_threshold).replace('.', '')+'_'+experiment+'_'+region2+'_'+model+'.png', bbox_inches='tight')


indxs3 = indxs3.astype(np.uint8)
indxs4 = indxs4.astype(np.uint8)
indxs5 = indxs5.astype(np.uint8)
indxs6 = indxs6.astype(np.uint8)

indxs3*= 1
indxs4*= 2
indxs5*= 3
indxs6*= 4

AS_TIP1 = np.ma.masked_where(indxs3==0, indxs3)
AS_TIP2 = np.ma.masked_where(indxs4==0, indxs4)
AS_TIP3 = np.ma.masked_where(indxs5==0, indxs5)
AS_TIP4 = np.ma.masked_where(indxs6==0, indxs6)

###################### save information in a mat dataset #####################
mdic = {'tips1': AS_TIP1, 'tips2': AS_TIP2, 'tips3': AS_TIP3, 'tips4': AS_TIP4}
savemat(path3+var+'_'+model+'_'+experiment+'_'+region2+'_tips_data.mat', mdic)


cmap8 = colors.ListedColormap([colours[0],colours[1],colours[2],colours[3]])

norm8 = colors.BoundaryNorm([0.5,1.5,2.5,3.5,4.5], cmap8.N)


# ##### Plot abrupt shift type #####

fig, ax = plt.subplots(figsize=figsizes)
fig, ax = plt.subplots(subplot_kw=dict(projection=proj),figsize=figsizes)
ax.set_extent([lonmin ,lonmax, latmin, latmax], crs=ccrs.PlateCarree())   
ax.add_feature(cfeature.COASTLINE, zorder=10)

ax.gridlines(draw_labels=False, ylocs = np.arange(-20,20,10), linestyle='--')

ax.set_xticks(np.arange(-80,-30,20),crs=ccrs.PlateCarree())
ax.set_xticklabels(np.arange(-80,-30,20))
ax.set_yticks(np.arange(-20,20,latstep),crs=ccrs.PlateCarree())
ax.set_yticklabels(np.arange(-20,20,latstep))
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.grid(linewidth=0.5, color='black', linestyle='--')
# ax.imshow(AS_TIP1 ,cmap=cmap8,norm=norm8,transform=ccrs.PlateCarree(),extent=[lonmin ,lonmax, latmax, latmin])
ax.imshow(AS_TIP2 ,cmap=cmap8,norm=norm8,transform=ccrs.PlateCarree(),extent=[lonmin ,lonmax, latmax, latmin])
# ax.imshow(AS_TIP3 ,cmap=cmap8,norm=norm8,transform=ccrs.PlateCarree(),extent=[lonmin ,lonmax, latmax, latmin])
# ax.imshow(AS_TIP4 ,cmap=cmap8,norm=norm8,transform=ccrs.PlateCarree(),extent=[lonmin ,lonmax, latmax, latmin])
for i in range (4):
    patches = [ mpatches.Patch(color=colours[i], label=type_labels[i] )]
fig.text(0.90, 0.71, 'a', ha='center')

# plt.legend(handles=patches, loc='lower center', bbox_to_anchor=(1, 0., 0.5, 0.5), frameon=False)  ## changed from loc=loc
# plt.legend(handles=patches, loc='upper center', bbox_to_anchor=(0.5, -0.05),
#           fancybox=True, shadow=True, ncol=5)
# plt.title(model)

plt.title (str(std_threshold)+' std deviations')
fig.tight_layout()

# if CONTROL:
#     fig.savefig(save_path+'/'+model+'_std_thresh'+str(std_threshold)+'.svg', format = 'svg', dpi=300, bbox_inches='tight')
#     # fig.savefig(save_path+'Control/'+var+'_as_type_asthresh'+str(abs_change_threshold)+'_fracthresh'+str(as_frac_threshold).replace('.', '')+'_stdthresh'+str(std_threshold)+'_'+experiment+'_'+region2+'_'+model+'_legend.svg', format = 'svg', dpi=300, bbox_inches='tight')
# else:
#     fig.savefig(save_path+var+'_as_type_asthresh'+str(abs_change_threshold)+'_fracthresh'+str(as_frac_threshold).replace('.', '')+'_'+experiment+'_'+region2+'_'+model+'_legend.png', bbox_inches='tight')

    
if CONTROL:
      ##### Plot std dev. multiple of control #####
    fig, ax = plt.subplots(figsize=figsizes)
    fig, ax = plt.subplots(subplot_kw=dict(projection=proj),figsize=figsizes)
    ax.set_extent([lonmin ,lonmax, latmin, latmax], crs=ccrs.PlateCarree())   
    ax.add_feature(cfeature.COASTLINE, zorder=10)
    ax.set_xticks(np.arange(-80,-30,lonstep),crs=ccrs.PlateCarree())
    ax.set_xticklabels(np.arange(-80,-30,lonstep))
    ax.set_yticks(np.arange(-20,20,latstep),crs=ccrs.PlateCarree())
    ax.set_yticklabels(np.arange(-20,20,latstep))
    lon_formatter = cticker.LongitudeFormatter()
    lat_formatter = cticker.LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.grid(linewidth=0.5, color='black', linestyle='--')
        
    im = ax.imshow(AS_STD,cmap=cmap2,norm=norm2,transform=ccrs.PlateCarree(),extent=[lonmin ,lonmax, latmax, latmin])
    cbar = fig.colorbar(im, pad=0.08, orientation='horizontal',extend='max')#, spacing="proportional")
    cbar.set_label('Ratio of abrupt shift gradient to std. dev. of control gradient')
    plt.title(model)
    # fig.tight_layout()
    # fig.savefig(save_path+'Control/'+var+'_as_std_asthresh'+str(abs_change_threshold)+'_fracthresh'+str(as_frac_threshold).replace('.', '')+'_stdthresh'+str(std_threshold)+'_'+experiment+'_'+region2+'_'+model+'.png', bbox_inches='tight')
