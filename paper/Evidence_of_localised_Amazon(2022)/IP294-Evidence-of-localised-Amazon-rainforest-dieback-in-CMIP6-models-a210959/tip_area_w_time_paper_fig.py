# -*- coding: utf-8 -*-
"""
Created on Thur Sept 30th 2021
Create figure for number of points tipped with increasing CO2

@author: impy2
"""
from pathlib import Path
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})
from scipy.io import loadmat
from scipy.io import savemat
import iris
import matplotlib.path as mpltPath

################# Specify the experiment and models ################

experiment = '1pctCO2'
variant_id = 'r1i1p1f1'

################# Specify variable of interest ######################
# var2 is the same as var but all lower case
var = 'cVeg' # or treeFrac
var2 = 'cveg'
var3 = 'tas'

############## Specify units of variable ###########################
units = '$kg/m^2$'

# Define latitude and longitude (uses 1 degree by 1 degree grid) 
my_lats=np.arange(-90,91,step=1)
my_lons=np.arange(-180,181,step=1)

ny =  np.size(my_lats)
nx = np.size(my_lons)

coord1 = iris.coords.DimCoord(my_lats,bounds=np.array([my_lats-0.5,my_lats+0.5]).T, standard_name='latitude', units='degrees', var_name='lat', attributes={'title': 'Latitude', 'type': 'double', 'valid_max': 90.0, 'valid_min': -90.0}, circular=True)

coord2 = iris.coords.DimCoord(my_lons,bounds=np.array([my_lons-0.5,my_lons+0.5]).T, standard_name='longitude', units='degrees', var_name='lon', attributes={'title': 'Longitude', 'type': 'double', 'valid_max': 180.0, 'valid_min': -180.0}, circular=True)

cube = iris.cube.Cube(np.zeros((ny,nx),np.float32),dim_coords_and_dims=[(coord1,0),(coord2,1)])

areas = iris.analysis.cartography.area_weights(cube, normalize=False)

## initialise figure
fig, axes = plt.subplots(1, 1, figsize=(2.7, 2.0), sharex=True)

################# Specify models here ############################# 
models = ['EC-Earth3-Veg', 'GFDL-ESM4', 'MPI-ESM1-2-LR', 'NorCPM1', 'SAM0-UNICON', 'TaiESM1', 'UKESM1-0-LL']
colors = ['tab:blue', 'tab:green', 'tab:purple', 'darkred', 'dimgray', 'goldenrod','lightseagreen']

## initialise zero arrays
comp_perc = np.zeros((np.size(models), 140))
comp_anom = np.zeros((np.size(models), 140))
wl = 15

# Is variable measured in per second?
PERSEC = False

## loop over each model
for v in range (0,np.size(models)):
    model = models[v]

    # model parameters
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
        variant_id = 'r1i1p1f2'
        
    ############ Specify name of directory of data which you want to plot #################
    region = 'Amazon'
    
    ############ Specify name of directory of that stores the figures #################
    region2 = 'Amazon'
    
    # Is variable measured in per second?
    PERSEC = False
        
    ######################## specify directory ###################################
    
    path = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment+'/Analysis_data/'+region2+'/'
    path2 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment+'/tip_points_data/'+region2+'/'
    path3 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment+'/Processed_data_monthly/'+region+'/'
    
    ####################### read in tip points data ##############################
    
    mat = loadmat(path+var+'_'+model+'_as_grads_wl'+str(wl)+'_data_'+experiment+'_'+region2+'.mat')
    mat2 = loadmat(path2+var+'_'+model+'_'+experiment+'_'+region+'_as_change_data.mat')
    
    ################ extract the time series for each tip point ##################
    # File name of interpolated data
    fname = path3+var+'_'+model+'_'+experiment+'_'+variant_id+'_'+date_range+'_mon_'+region+'.nc'
    fname = Path(fname)
    
    # # Load in data
    f = nc.Dataset(fname,'r')
    # Extract data and if measured in per second convert to per year
    if PERSEC:
        x = f.variables[var2][:]*3600*24*360
    else:
        x = f.variables[var2][:]
    #close dataset
    f.close()
    # Extract data from dictionary
    Lon, Lat, as_grad, as_grads, co2_tip, as_change, ovr_change = mat['Lon'], mat['Lat'], mat['as_grad'], mat['as_grads'], mat['co2_tip'], mat['as_change'], mat['ovr_change']
         
    ### define tipping point array ####
    tips = mat2['tips']
    
    #################### create masks of the Amazon regions ######################
    lon, lat = np.meshgrid(Lon[0,:], Lat[:,0])
    Lon_flat, Lat_flat = Lon.flatten(), Lat.flatten()
    points = np.vstack((Lon_flat,Lat_flat)).T
    
    # Vertices of the 4 Amazon boxes
    # NWS_verts = [(-75,12), (-83.4,2.2), (-83.4,-10), (-79,-15), (-72,-15), (-72,12)]
    NSA_verts = [(-72,12), (-72,-8), (-50,-8), (-50,7.6), (-55,12)]
    # NES_verts = [(-34,-20), (-50,-20), (-50,0), (-34,0)]
    # SAM_verts = [(-66.4,-20), (-72,-15), (-72,-8), (-50,-8), (-50,-20)]
    
    # nt = np.size(tips, 0)
    ny = np.size(tips, 0)
    nx = np.size(tips, 1)
    
    # Create masks for each region (note nt, ny, nx, are dimension sizes of time, latitude, longitude respectively)
    # NWS_path = mpltPath.Path(NWS_verts)
    # NWS_grid = NWS_path.contains_points(points)
    # NWS_grid = NWS_grid.reshape((ny, nx))
    # NWS_mask = np.broadcast_to(NWS_grid, (ny, nx))
    
    NSA_path = mpltPath.Path(NSA_verts)
    NSA_grid = NSA_path.contains_points(points)
    NSA_grid = NSA_grid.reshape((ny, nx))
    NSA_mask = np.broadcast_to(NSA_grid, (ny, nx))
    
    # NES_path = mpltPath.Path(NES_verts)
    # NES_grid = NES_path.contains_points(points)
    # NES_grid = NES_grid.reshape((ny, nx))
    # NES_mask = np.broadcast_to(NES_grid, (ny, nx))
    
    # SAM_path = mpltPath.Path(SAM_verts)
    # SAM_grid = SAM_path.contains_points(points)
    # SAM_grid = SAM_grid.reshape((ny, nx))
    # SAM_mask = np.broadcast_to(SAM_grid, (ny, nx))
    
    # Create copies of surface temperature for each Amazon region
    # NWS_tips = tips.copy()
    NSA_tips = tips.copy()
    # NES_tips = tips.copy()
    # SAM_tips = tips.copy()
    
    NSA_lat = Lat.astype(float)
    
    # Mask out all data not in Amazon box
    # NWS_tips[~NWS_mask] = np.nan
    NSA_tips[~NSA_mask] = np.nan
    # NES_tips[~NES_mask] = np.nan
    # SAM_tips[~SAM_mask] = np.nan
    
    NSA_lat[~NSA_mask] = np.nan
    
    ######### extract indicies of tipping points ########################
    [rows, cols] = np.where(NSA_tips<1e+20)
    
    ######## calcualte the areas corresponding to these grid points ########
    ## the radius of the Earth, r ##
    r = 6357 # (km)
    
    ## define a function to calculate the area of each grid point
    def f_area(mid_lat):
        mid_lat_r = mid_lat*(np.pi/180)
        A = (((2* np.pi *r)/360)**2)*(np.cos(mid_lat_r))
        return A
    
    grid_areas = np.zeros(np.size(rows))
    for i in range (0, np.size(rows)):
        latitude = Lat[rows[i], cols[i]]
        grid_areas[i] = f_area(latitude) 
    
    am_areas = np.zeros((34, 49))
    for i in range(0, 34):
        for j in range(0,49):
            latitude = NSA_lat[i,j]
            am_areas[i,j] = f_area(latitude)
    
    #### total area of amazon
    am_areas = am_areas.flatten()
    am_areas = am_areas[~np.isnan(am_areas)]
    total_a = np.sum(am_areas)
    
    ###################### extract tipping points from data #######################
    
    co2_start = 284.3186666723341
    co2 = np.zeros(total_years)
    
    for k in range(total_years):
        co2[k] = co2_start*(1.01**k)
        
    co2_tip = mat['co2_tip']
    
    #### create array of co2 at each tipping point #####
    co2_at_tip = np.zeros(np.size(rows))
    for i in range (0, np.size(rows)):
        co2_at_tip[i] = co2_tip[rows[i], cols[i]]
    
    count = np.zeros(np.size(co2))
    
    ####### calcualte the number of tip points that have tipped at each co2 level ######
    for i in range (0, np.size(co2)):
        level = co2[i]
        for k in range (0, np.size(rows)):
            area = grid_areas[k]
            tip_level = co2_at_tip[k]
            if tip_level <= level:
                count[i] = count[i] + area
    
    #### calculate the percentage ####            
    perc = np.zeros(np.size(count))            
    for i in range (np.size(count)):
        perc[i] = count[i]/total_a 
    
    
    ############### calcualte the temperature anomolies over time ################
    path4 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var3+'/'+experiment+'/Processed_data_monthly/'
    region1 = 'World'
    
    # File name of interpolated World data
    fname = path4+region1+'/'+var3+'_'+model+'_'+experiment+'_'+variant_id+'_'+date_range+'_mon_'+region1+'.nc'
    # Load in world data into iris cube
    world_cube = iris.load_cube(fname,var3)
    
    world_data = world_cube.data
    
    #################### average temperature across the globe ####################
    
    weighted_tas =np.zeros((np.size(world_data,0), np.size(world_data,1), np.size(world_data,2)))
    weighting = areas/np.max(areas)
    for i in range (0, np.size(world_data,0)):
        weighted_tas[i] = world_data[i,:,:]*weighting
        
    pavr = np.zeros(np.size(weighted_tas, 0))
    for i in range (np.size(weighted_tas, 0)):
        t_point = weighted_tas[i,:,:]
        t_point = t_point.flatten()
        t_point = t_point[~np.isnan(t_point)]
        pavr[i] = np.average(t_point)
    
        tavr = np.zeros(int(np.size(weighted_tas, 0)/12))
    for i in range (np.size(tavr)):
        tavr[i] = np.average(pavr[i*12:(i*12)+12])          
        
    ### first 10 years of the model is the ref temp
    ref = np.average(tavr[0:10])
    
    anom = np.zeros(np.size(tavr)-10)
    #### calculate the temperaure anomaly for each 10 year window going forward
    for i in range (0, np.size(tavr)-10):
        temp = np.average(tavr[i:i+10])
        anom[i] = temp - ref
    
    comp_perc[v, :] = perc[10:150]
    comp_anom[v,:] = anom[0:140]
    
    ### plot model line on graph ###
    axes.plot(anom, perc[10::]*100, label=model, color=colors[v], linewidth=0.8)

################## save data as a mat file ################################

# Create dictionary of data to save                
mdic = {'EC_Earth3-Veg':comp_perc[0,:], 'GFDL_ESM4':comp_perc[1,:], 'MPI-ESM1-2-LR': comp_perc[2,:], 'NorCPM1': comp_perc[3,:], 'SAM0-UNCON': comp_perc[4,:], 'TaiESM1': comp_perc[5,:], 'UKESM1-0-LL':comp_perc[6,:]}
mdic2 = {'EC_Earth3-Veg':comp_anom[0,:], 'GFDL_ESM4':comp_anom[1,:], 'MPI-ESM1-2-LR': comp_anom[2,:], 'NorCPM1': comp_anom[3,:], 'SAM0-UNCON': comp_anom[4,:], 'TaiESM1': comp_anom[5,:], 'UKESM1-0-LL':comp_anom[6,:]}
savepath_mat = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/Figures/Paper/Figure 2/'

## Save data as a mat fifle - Specify filename ##
savemat(savepath_mat+'shift_area_w_warming', mdic)
savemat(savepath_mat+'global_anoms', mdic2)


#### create x axis coordinates to interpolate onto 
x_interp = np.arange(0, 4, step = 4/np.size(comp_anom,1))

#### interpolate the data ####
interp_perc = np.zeros((np.size(comp_perc,0), np.size(comp_perc, 1)))
for i in range (0, np.size(interp_perc,0)):
    interp_perc[i,:] = np.interp(x_interp, comp_anom[i,:], comp_perc[i,:])

compiled_perc = np.zeros(140)
for i in range (0, np.size(compiled_perc)):
    compiled_perc[i] = np.sum(interp_perc[:, i])/np.size(models)
    
axes.plot(x_interp, compiled_perc*100, label='Compiled Models', color='black', linewidth=2.0)

##### calculate error in each index of the array ######
err = np.zeros(np.size(interp_perc,1))
for i in range (0,np.size(interp_perc,1)):
    err[i] = np.std(interp_perc[:,i]*100)

up_err = (compiled_perc*100)+err
low_err = (compiled_perc*100)-err

lower_bound = np.zeros(np.size(low_err))
for i in range (0, np.size(low_err)):
    lower_bound[i]= np.maximum(0, low_err[i])

axes.fill_between(x_interp, lower_bound, up_err, alpha=0.25, color='black')
        
plt.ylabel ('Area with Abrupt \n Shifts (%)')
plt.xlabel('Global Warming ('+ u'\u2103'+')')
plt.xticks(np.arange(0, 4.5, 0.5))
plt.xlim([0,3.0])
plt.yticks(np.arange(0, 50, 10))
fig.text(0.85, 0.79, '(b)', ha='center')
plt.show()

# savepath = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/Figures/Paper/Figure 2/'
# filename = savepath+'area_tipped_V7.svg'
# fig.savefig(filename, format = 'svg', dpi=300, bbox_inches='tight') 

    
    
    
    