# -*- coding: utf-8 -*-
"""
Created on Thu May 20 16:54:16 2021

Boostrap of linear regression analysis found in Fig 4

@author: impy2
"""

from pathlib import Path
import netCDF4 as nc
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 10})

from scipy.io import loadmat
import iris
from scipy.io import savemat
import matplotlib.path as mpltPath
################# Specify the experiment and models ################

experiment = '1pctCO2'
experiment2 = 'piControl'


################# Specify variable of interest ######################
var = 'tas' # or treeFrac
var2 = 'tas' # or treefrac
tip_var = 'cVeg'
tip_var2 = 'cveg'

############## Specify units of variable ###########################
units = '$kgm^{-2}$'

################ Specify model here ################################## 
models = ['EC-Earth3-Veg', 'GFDL-ESM4', 'MPI-ESM1-2-LR', 'NorCPM1', 'SAM0-UNICON', 'TaiESM1', 'UKESM1-0-LL']
# models = ['NorCPM1']
colors = ['tab:blue', 'tab:green', 'tab:purple', 'darkred', 'dimgray', 'goldenrod', 'lightseagreen']
alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']

#define bins for histograms
thresholds = np.arange(-0.6, 3.1, 0.2)
risk = np.zeros((np.size(models),np.size(thresholds)))

# Define latitude and longitude (uses 1 degree by 1 degree grid) 
my_lats=np.arange(-90,91,step=1)
my_lons=np.arange(-180,181,step=1)

ny =  np.size(my_lats)
nx = np.size(my_lons)

coord1 = iris.coords.DimCoord(my_lats,bounds=np.array([my_lats-0.5,my_lats+0.5]).T, standard_name='latitude', units='degrees', var_name='lat', attributes={'title': 'Latitude', 'type': 'double', 'valid_max': 90.0, 'valid_min': -90.0}, circular=True)

coord2 = iris.coords.DimCoord(my_lons,bounds=np.array([my_lons-0.5,my_lons+0.5]).T, standard_name='longitude', units='degrees', var_name='lon', attributes={'title': 'Longitude', 'type': 'double', 'valid_max': 180.0, 'valid_min': -180.0}, circular=True)

cube = iris.cube.Cube(np.zeros((ny,nx),np.float32),dim_coords_and_dims=[(coord1,0),(coord2,1)])

areas = iris.analysis.cartography.area_weights(cube, normalize=False)

# where 'point' is the time at which 73 years have passed and CO2 has doubled
point = 876
time = np.arange(0, 73, 1)

# Create array of CO2 for 1pctCO2 run
co2_start = 284.3186666723341
co2 = np.zeros(int(point/12))
for k in range(int(point/12)):
    co2[k] = co2_start*(1.01**k)
    
log_co2 = np.log10(co2)    

## define a function to calculate the area of each grid point
## the radius of the Earth, r ##
r = 6357 # (km)

def f_area(mid_lat):
    mid_lat_r = mid_lat*(np.pi/180)
    A = (((2* np.pi *r)/360)**2)*(np.cos(mid_lat_r))
    return A
   
#start figure 
fig, ax = plt.subplots(3, 3, figsize=(8.27, 6.25))
plt.tight_layout()
ax = ax.flatten()

#%%

##loop for each model removed
for q in range (0,np.size(models)):
# for q in range(0,1):
    
    #define model lists
    models = ['EC-Earth3-Veg', 'GFDL-ESM4', 'MPI-ESM1-2-LR', 'NorCPM1', 'SAM0-UNICON', 'TaiESM1', 'UKESM1-0-LL']
    models1 = ['EC-Earth3-Veg', 'GFDL-ESM4', 'MPI-ESM1-2-LR', 'NorCPM1', 'SAM0-UNICON', 'TaiESM1', 'UKESM1-0-LL']
    #define legend labels 
    models2 = ['- EC-Earth3-Veg', '- GFDL-ESM4', '- MPI-ESM1-2-LR', '- NorCPM1', '- SAM0-UNICON', '- TaiESM1', '- UKESM1-0-LL']

    models.remove(models[q])
    ######### define arrays in which final inforation will be stored  ############
    tip1 = np.zeros((np.size(models), 1666))
    tip2 = np.zeros((np.size(models), 1666))
    tip3 = np.zeros((np.size(models), 1666))
    tip4 = np.zeros((np.size(models), 1666))
    notip = np.zeros((np.size(models), 1666))
    
    tot_area1 = np.zeros((np.size(models), 1666))
    tot_area2 = np.zeros((np.size(models), 1666))
    tot_area3 = np.zeros((np.size(models), 1666))
    tot_area4 = np.zeros((np.size(models), 1666))
    tot_arean = np.zeros((np.size(models), 1666))
    
    ###perform analysis with one model missing ####
    for l in range (0, np.size(models)):
        
        model = models[l]
        
        if model == 'UKESM1-0-LL':
            variant_id = 'r1i1p1f2'
        else:
            variant_id = 'r1i1p1f1'
        # Specific entries for the 3 models initially used
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
        
        ############### calcualte the temperature anomolies over time ################
        path4 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var2+'/'+experiment+'/Processed_data_monthly/'
        region1 = 'World'
        
        # File name of interpolated World data
        fname = path4+region1+'/'+var2+'_'+model+'_'+experiment+'_'+variant_id+'_'+date_range+'_mon_'+region1+'.nc'
        # Load in world data into iris cube
        world_cube = iris.load_cube(fname,var2)
        
        world_data = world_cube.data
        
        #################### average temperature across the globe ####################
        
        weighted_tas =np.zeros((np.size(world_data,0), np.size(world_data,1), np.size(world_data,2)))
        weighting = areas/np.max(areas)
        for i in range (0, np.size(world_data,0)):
            weighted_tas[i] = world_data[i,:,:]*weighting
        
        #################### average temperature across the globe ####################
        
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
        
        x_array = anom[0:73]    
        
        
        # Is variable measured in per second?
        PERSEC = False
        
        CONTROL = True
            
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
        path5 =  'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+'cVeg'+'/'+experiment+'/Figure_data/'
        
        if CONTROL:
            ################## Specify path to analysis data directory #####################
            path4 = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+tip_var+'/'+experiment2+'/Analysis_data/'+region2+'/'
        
            ################## Specify filename of saved analysis data and load in #####################
            mat2 = loadmat(path4+tip_var+'_'+model+'_wl'+str(wl)+'_data_'+experiment2+'_'+region2+'.mat')
        
            control_grads = mat2['control_grads']
            control_grad_std = np.nanstd(control_grads[-400:-wl,:,:], axis=0)
        
            ############ Set minimum number of standard deviations from ############
            ############ zero for it to be classed as an abrupt shift ############
            std_threshold = 3
        
        # File name of interpolated data
        fname = path+var+'_'+model+'_'+experiment+'_'+variant_id+'_'+date_range+'_mon_'+region+'.nc'
        fname = Path(fname)
        
        # # Load in data
        f = nc.Dataset(fname,'r')
        
        mat = loadmat(path2+tip_var+'_'+model+'_as_grads_wl'+str(wl)+'_data_'+experiment+'_'+region2+'.mat')
        
        # Extract data from dictionary
        longitude, latitude, as_grad, as_grads, co2_tip, as_change, ovr_change = mat['Lon'], mat['Lat'], mat['as_grad'], mat['as_grads'], mat['co2_tip'], mat['as_change'], mat['ovr_change']
        
        ############ Set minimum threshold for change in sliding window to #############
        ############ be classed as an abrupt shift #####################################
        abs_change_threshold = 2
        
        ############ Set minimum fractional contribution of abrupt shift to ##########
        ############ overall change ##################################################
        as_frac_threshold = 0.25
        
        # Extract data and if measured in per second convert to per year
        if PERSEC:
            x = f.variables[var2][:]*3600*24*360
        else:
            x = f.variables[var2][:]
        
        # Close dataset
        f.close()
        
        #################### create masks of the Amazon regions ######################
        Lon, Lat = np.meshgrid(longitude[0,:], latitude[:,0])
        Lon_flat, Lat_flat = Lon.flatten(), Lat.flatten()
        points = np.vstack((Lon_flat,Lat_flat)).T
        
        # Vertices of the 4 Amazon boxes
        NWS_verts = [(-75,12), (-83.4,2.2), (-83.4,-10), (-79,-15), (-72,-15), (-72,12)]
        NSA_verts = [(-72,12), (-72,-8), (-50,-8), (-50,7.6), (-55,12)]
        NES_verts = [(-34,-20), (-50,-20), (-50,0), (-34,0)]
        SAM_verts = [(-66.4,-20), (-72,-15), (-72,-8), (-50,-8), (-50,-20)]
        
        nt = np.size(x, 0)
        ny = np.size(x, 1)
        nx = np.size(x, 2)
        
        # Create masks for each region (note nt, ny, nx, are dimension sizes of time, latitude, longitude respectively)
        NWS_path = mpltPath.Path(NWS_verts)
        NWS_grid = NWS_path.contains_points(points)
        NWS_grid = NWS_grid.reshape((ny, nx))
        NWS_mask = np.broadcast_to(NWS_grid, (nt, ny, nx))
        
        NSA_path = mpltPath.Path(NSA_verts)
        NSA_grid = NSA_path.contains_points(points)
        NSA_grid = NSA_grid.reshape((ny, nx))
        NSA_mask = np.broadcast_to(NSA_grid, (nt, ny, nx))
        
        NES_path = mpltPath.Path(NES_verts)
        NES_grid = NES_path.contains_points(points)
        NES_grid = NES_grid.reshape((ny, nx))
        NES_mask = np.broadcast_to(NES_grid, (nt, ny, nx))
        
        SAM_path = mpltPath.Path(SAM_verts)
        SAM_grid = SAM_path.contains_points(points)
        SAM_grid = SAM_grid.reshape((ny, nx))
        SAM_mask = np.broadcast_to(SAM_grid, (nt, ny, nx))
        
        # Create copies of surface temperature for each Amazon region
        NWS_tas = x.copy()
        NSA_tas = x.copy()
        NES_tas = x.copy()
        SAM_tas = x.copy()
        
        # Mask out all data not in Amazon box
        NWS_tas[~NWS_mask] = np.nan
        NSA_tas[~NSA_mask] = np.nan
        NES_tas[~NES_mask] = np.nan
        SAM_tas[~SAM_mask] = np.nan
        
        #define 4 different types of tipping point
        indxs = np.where((np.abs(as_grad)>std_threshold*control_grad_std)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
        indxs2 = ~((np.abs(as_grad)>std_threshold*control_grad_std)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
        indxs3 = ((ovr_change>0)&(as_change>0)&(np.abs(as_grad)>std_threshold*control_grad_std)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
        indxs4 = ((ovr_change<0)&(as_change<0)&(np.abs(as_grad)>std_threshold*control_grad_std)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
        indxs5 = ((ovr_change<0)&(as_change>0)&(np.abs(as_grad)>std_threshold*control_grad_std)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
        indxs6 = ((ovr_change>0)&(as_change<0)&(np.abs(as_grad)>std_threshold*control_grad_std)&(np.abs(as_change)>abs_change_threshold)&(np.abs(as_change/ovr_change)>as_frac_threshold))
        
            
        indxs3 = indxs3.astype(np.uint8)
        indxs4 = indxs4.astype(np.uint8)
        indxs5 = indxs5.astype(np.uint8)
        indxs6 = indxs6.astype(np.uint8)
        
        indxs3*= 1
        indxs4*= 2
        indxs5*= 3
        indxs6*= 4
        
        #create masked arrays for each tipping point 
        AS_TIP1 = np.ma.masked_where(indxs3==0, indxs3)
        AS_TIP2 = np.ma.masked_where(indxs4==0, indxs4)
        AS_TIP3 = np.ma.masked_where(indxs5==0, indxs5)
        AS_TIP4 = np.ma.masked_where(indxs6==0, indxs6)
        
        #masked array for all tipping points 
        AS_CHANGE = np.ma.masked_where(indxs2, as_change)
                
        year_stop = 73
        
        #create zero arrays for seasonal cycle amplitudes 
        diff_y1 = np.zeros((year_stop, np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        diff_y2 = np.zeros((year_stop, np.size(AS_TIP2,0), np.size(AS_TIP2,1)))
        diff_y3 = np.zeros((year_stop, np.size(AS_TIP3,0), np.size(AS_TIP3,1)))
        diff_y4 = np.zeros((year_stop, np.size(AS_TIP4,0), np.size(AS_TIP4,1)))
        diff_yn = np.zeros((year_stop, np.size(AS_CHANGE,0), np.size(AS_CHANGE,1)))
        
        #create zero arrays for ratios 
        # ratio1 = np.zeros((year_stop, np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        # ratio2 = np.zeros((year_stop, np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        # ratio3 = np.zeros((year_stop, np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        # ratio4 = np.zeros((year_stop, np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        # ration = np.zeros((year_stop, np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        
        #create zero arrays for linear regression outputs 
        b0_1 = np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        a0_1 = np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        b0_2 = np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        a0_2 = np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        b0_3 = np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        a0_3 = np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        b0_4 = np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        a0_4 = np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        b0_n = np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        a0_n = np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        
        #### arrays for the areas of each var ###
        area1 =np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1))) 
        area2 =np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1))) 
        area3 =np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1))) 
        area4 =np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        arean =np.zeros((np.size(AS_TIP1,0), np.size(AS_TIP1,1)))
        
        ########### calculate kendall tau values for each tip point type #############
        #### tip point type 1 ####
        for i in range (0, np.size(AS_TIP1,0)):
            for j in range (0,np.size(AS_TIP1,1)):
                d = AS_TIP1[i,j]
                if d == 1:
                    T = NSA_tas[:,i,j].flatten()
                    if np.isnan(T[0]) == False:                    
                        T = T[0:point]
                        #### separate data into years #####
                        s = int(np.size(T)/12)
                        year_points = np.zeros((s, 12))
                        max_y = np.zeros(s)
                        min_y = np.zeros(s)
                                   
                        for k in range (0, s):
                            year_points[k, :] = np.arange(k*12, (k*12)+12)
                        
                        #extract data per year and then find the max temp in each yr 
                        for k in range (0,s):
                            points = year_points[k,:]
                            d_y = T[int(points[0]):int(points[-1])+1]
                            max_y[k] = np.max(d_y)
                            min_y[k] = np.min(d_y)
                            diff_y1[k,i,j] = max_y[k] - min_y[k]
                            # ratio1[k,i,j] = (diff_y1[k,i,j]-diff_y1[0,i,j])/diff_y1[0,i,j]
                        
                        ##### calculate the linear regression for this point ####
                        b0_1[i,j], a0_1[i,j], r0, p0, e0 = stats.linregress(x_array,diff_y1[:,i,j])
                        # b0_1[i,j], a0_1[i,j], r0, p0, e0 = stats.linregress(x_array,ratio1[:,i,j])
                        
                        ### calculate corresponding area ###
                        a_lat = latitude[i,j]
                        area1[i,j] = f_area(a_lat)
                    else:
                        area1[i,j] = np.nan
                        b0_1[i,j] = np.nan
        
                
        b0_1 = b0_1.flatten()
        tip1[l,:] = b0_1
        b0_1 = b0_1[~np.isnan(b0_1)]
        b0_1 = b0_1[b0_1 != 0]
        
        area1 = area1.flatten()
        tot_area1[l,:] = area1
        area1 = area1[~np.isnan(area1)]
        area1 = area1[area1 != 0]
    
        
        #### for tip point 2 ####
        for i in range (0, np.size(AS_TIP2,0)):
            for j in range (0,np.size(AS_TIP2,1)):
                d = AS_TIP2[i,j]
                if d == 2:
                    T = NSA_tas[:,i,j].flatten()
                    if np.isnan(T[0]) == False:
                        T = T[0:point]
                        #### separate data into years #####
                        s = int(np.size(T)/12)
                        year_points = np.zeros((s, 12))
                        max_y = np.zeros(s)
                        min_y = np.zeros(s)
                                   
                        for k in range (0, s):
                            year_points[k, :] = np.arange(k*12, (k*12)+12)
                        
                        #extract data per year and then find the max temp in each yr 
                        for k in range (0,s):
                            points = year_points[k,:]
                            d_y = T[int(points[0]):int(points[-1])+1]
                            max_y[k] = np.max(d_y)
                            min_y[k] = np.min(d_y)
                            diff_y2[k,i,j] = max_y[k] - min_y[k] 
                            # ratio2[k,i,j] = diff_y2[k,i,j]/diff_y2[0,i,j]
                            
                        ##### calculate the linear regression for this point ####
                        b0_2[i,j], a0_2[i,j], r0, p0, e0 = stats.linregress(x_array,diff_y2[:,i,j])
                        # b0_2[i,j], a0_2[i,j], r0, p0, e0 = stats.linregress(x_array,ratio2[:,i,j])
                        
                        ### calculate corresponding area ###
                        a_lat = latitude[i,j]
                        area2[i,j] = f_area(a_lat)
                    else:
                        area2[i,j] = np.nan
                        b0_2[i,j] = np.nan
        
        where_b02 = b0_2            
        b0_2 = b0_2.flatten()
        tip2[l,:] = b0_2
        b0_2 = b0_2[~np.isnan(b0_2)]
        b0_2 = b0_2[b0_2 != 0]
    
        area2 = area2.flatten()
        tot_area2[l,:] = area2
        area2 = area2[area2 != 0]
        area2 = area2[~np.isnan(area2)]

        #### for tip point 3 ####
        for i in range (0, np.size(AS_TIP3,0)):
            for j in range (0,np.size(AS_TIP3,1)):
                d = AS_TIP3[i,j]
                if d == 3:
                    T = NSA_tas[:,i,j].flatten()
                    if np.isnan(T[0]) == False:        
                        T = T[0:point]
                        #### separate data into years #####
                        s = int(np.size(T)/12)
                        year_points = np.zeros((s, 12))
                        max_y = np.zeros(s)
                        min_y = np.zeros(s)
                                   
                        for k in range (0, s):
                            year_points[k, :] = np.arange(k*12, (k*12)+12)
                        
                        #extract data per year and then find the max temp in each yr 
                        for k in range (0,s):
                            points = year_points[k,:]
                            d_y = T[int(points[0]):int(points[-1])+1]
                            max_y[k] = np.max(d_y)
                            min_y[k] = np.min(d_y)
                            diff_y3[k,i,j] = max_y[k] - min_y[k]
                            # ratio3[k,i,j] = (diff_y3[k,i,j]-diff_y3[0,i,j])/diff_y3[0,i,j]
                            
                        ##### calculate the linear regression for this point ####
                        b0_3[i,j], a0_3[i,j], r0, p0, e0 = stats.linregress(x_array,diff_y3[:,i,j])
                        # b0_3[i,j], a0_3[i,j], r0, p0, e0 = stats.linregress(x_array,ratio3[:,i,j])
                        
                        ### calculate corresponding area ###
                        a_lat = latitude[i,j]
                        area3[i,j] = f_area(a_lat)
                    else:
                        b0_3[i,j] = np.nan
                        area3[i,j] = np.nan
                    
        b0_3 = b0_3.flatten()
        tip3[l,:] = b0_3
        b0_3 = b0_3[~np.isnan(b0_3)]
        b0_3 = b0_3[b0_3 != 0]
    
        area3 = area3.flatten()
        tot_area3[l,:] = area3
        area3 = area3[area3 != 0]
        
        #### for tip point 4 ####
        for i in range (0, np.size(AS_TIP4,0)):
            for j in range (0,np.size(AS_TIP4,1)):
                d = AS_TIP4[i,j]
                if d == 4:
                    T = NSA_tas[:,i,j].flatten()
                    if np.isnan(T[0]) == False:
                        T = T[0:point]
                        #### separate data into years #####
                        s = int(np.size(T)/12)
                        year_points = np.zeros((s, 12))
                        max_y = np.zeros(s)
                        min_y = np.zeros(s)
                                   
                        for k in range (0, s):
                            year_points[k, :] = np.arange(k*12, (k*12)+12)
                        
                        #extract data per year and then find the max temp in each yr 
                        for k in range (0,s):
                            points = year_points[k,:]
                            d_y = T[int(points[0]):int(points[-1])+1]
                            max_y[k] = np.max(d_y)
                            min_y[k] = np.min(d_y)
                            diff_y4[k,i,j] = max_y[k] - min_y[k]
                            # ratio4[k,i,j] = diff_y4[k,i,j]/diff_y4[0,i,j]
                            
                        ##### calculate the linear regression for this point ####
                        b0_4[i,j], a0_4[i,j], r0, p0, e0 = stats.linregress(x_array,diff_y4[:,i,j])
                        # b0_4[i,j], a0_4[i,j], r0, p0, e0 = stats.linregress(x_array,ratio4[:,i,j])
                         
                        ### calculate corresponding area ###
                        a_lat = latitude[i,j]
                        area4[i,j] = f_area(a_lat)
                    else:
                        area4[i,j] = np.nan
                        b0_4[i,j] = np.nan
                    
        b0_4 = b0_4.flatten()
        tip4[l,:] = b0_4
        b0_4 = b0_4[~np.isnan(b0_4)]
        b0_4 = b0_4[b0_4 != 0]
    
        area4 = area4.flatten()
        tot_area4[l,:] = area4
        area4 = area4[area4 != 0]
        area4 = area4[~np.isnan(area4)]
        
        ############ for the non tipping points ##############
        
        for i in range (0, np.size(AS_CHANGE,0)):
            for j in range (0, np.size(AS_CHANGE,1)):
                #identify if point is one of the four types of tipping point
                #AS_TIP masks are True where there is no tipping point detected 
                b = AS_TIP1[i,j]
                if np.ma.is_masked(b) is True:
                    b = AS_TIP2[i,j]
                    if np.ma.is_masked(b) is True:
                        b = AS_TIP3[i,j]
                        if np.ma.is_masked(b) is True:
                            b= AS_TIP4[i,j]
                            if np.ma.is_masked(b) is True:
                                if np.isfinite(as_grad[i,j]):
                                    T = NSA_tas[:,i,j].flatten()
                                    if np.isnan(T[0]) == False:
                                        T = T[0:point]
                                        #### separate data into years #####
                                        s = int(np.size(T)/12)
                                        year_points = np.zeros((s, 12))
                                        max_y = np.zeros(s)
                                        mean_y = np.zeros(s)
                                                   
                                        for k in range (0, s):
                                            year_points[k, :] = np.arange(k*12, (k*12)+12)
                                        
                                        #extract data per year and then find the max temp in each yr 
                                        for k in range (0,s):
                                            points = year_points[k,:]
                                            d_y = T[int(points[0]):int(points[-1])+1]
                                            max_y[k] = np.max(d_y)
                                            mean_y[k] = np.average(d_y)
                                            diff_yn[k,i,j] = max_y[k] - mean_y[k]
                                            # ration[k,i,j] = diff_yn[k,i,j]/diff_yn[0,i,j]
                                        
                                        ##### calculate the linear regression for this point ####
                                        b0_n[i,j], a0_n[i,j], r0, p0, e0 = stats.linregress(x_array,diff_yn[:,i,j])
                                        # b0_n[i,j], a0_n[i,j], r0, p0, e0 = stats.linregress(x_array,ration[:,i,j])
                                        
                                        ### calculate corresponding area ###
                                        a_lat = latitude[i,j]
                                        arean[i,j] = f_area(a_lat)
                                    else:
                                        b0_n[i,j] = np.nan
                                        arean[i,j] = np.nan 
        
        b0_n = b0_n.flatten()
        notip[l,:] = b0_n
        b0_n = b0_n[~np.isnan(b0_n)]
        b0_n = b0_n[b0_n != 0]    
    
        arean = arean.flatten()
        tot_arean[l,:] = arean
        arean = arean[arean != 0]   
        arean = arean[~np.isnan(arean)]                         
    
        ###### calcualte the risk of a tip point above a certain threshold ########
        count_tip = np.zeros(np.size(thresholds))
        count_n = np.zeros(np.size(thresholds))
        
        for i in range (0, np.size(thresholds)):
            thresh = thresholds[i]
            
            tipped = np.where(b0_2 < thresh)
            tipped_area = area2[tipped]
            
            count_tip[i] = np.sum(tipped_area)
    
        am_area = np.sum(area1) + np.sum(area2) + np.sum(area3) + np.sum(area4) + np.sum(arean)
    
    #%%   
    ########## create a line for the compliled models #############
    tot_area1 = tot_area1.flatten()
    tot_area1 = tot_area1[~np.isnan(tot_area1)]
    tot_area1 = tot_area1[tot_area1 !=0]
    
    tot_area2 = tot_area2.flatten()
    tot_area2 = tot_area2[~np.isnan(tot_area2)]
    tot_area2 = tot_area2[tot_area2 !=0]
    
    tot_area3 = tot_area3.flatten()
    tot_area3 = tot_area3[~np.isnan(tot_area3)]
    tot_area3 = tot_area3[tot_area3 !=0]
    
    tot_area4 = tot_area4.flatten()
    tot_area4 = tot_area4[~np.isnan(tot_area4)]
    tot_area4 = tot_area4[tot_area4 !=0]
    
    tot_arean = tot_arean.flatten()
    tot_arean = tot_arean[~np.isnan(tot_arean)]
    tot_arean = tot_arean[tot_arean !=0]
    
    areafrac1 = tot_area1/am_area
    areafrac2 = tot_area2/am_area
    areafrac3 = tot_area3/am_area
    areafrac4 = tot_area4/am_area
    areafracn = tot_arean/am_area
    
    ##### do the same for area #####
    tip1 = tip1.flatten()
    tip1 = tip1[~np.isnan(tip1)]
    tip1 = tip1[tip1 !=0]
    
    tip2 = tip2.flatten()
    tip2 = tip2[~np.isnan(tip2)]
    tip2 = tip2[tip2 !=0]
    
    tip3 = tip3.flatten()
    tip3 = tip3[~np.isnan(tip3)]
    tip3 = tip3[tip3 !=0]
    
    tip4 = tip4.flatten()
    tip4 = tip4[~np.isnan(tip4)]
    tip4 = tip4[tip4 !=0]
    
    notip = notip.flatten()
    notip = notip[~np.isnan(notip)]
    notip = notip[notip !=0]
    
    
    ########################### plot histogram #############################
    edges = thresholds
    s = np.size(models)
    ax[q].title.set_text('minus '+ models1[q])
    ax[q].hist([tip2, tip1, tip3, tip4, notip], weights = [(areafrac2/s)*100, (areafrac1/s)*100, (areafrac3/s)*100, (areafrac4/s)*100, (areafracn/s)*100], 
               stacked=True,color=['red','purple', 'purple','purple', 'purple'], bins = edges)
    ax[q].text(0.88, 0.88, '('+alphabet[q]+')', transform=ax[q].transAxes)
    ax[q].set_xticks(np.arange(-0.5, 3.5, 1))
#%%
    tip_count = np.zeros(np.size(thresholds))
    
    #### caclulate bar chart of risk for panel 8
    for i in range (1, np.size(thresholds)):
        low_thresh = thresholds[i-1]
        up_thresh = thresholds[i]
        ## for red coded tipping points 
        wheretip2 = np.logical_and(tip2 >= low_thresh, tip2 <= up_thresh)
        tipped = np.where(wheretip2 == True)
        tipped_area = np.sum(tot_area2[tipped])
        
        ## all other points 
        wheretip1 =np.logical_and(tip1 >= low_thresh, tip1 <= up_thresh)
        wheretip3 =np.logical_and(tip3 >= low_thresh, tip3 <= up_thresh)
        wheretip4 =np.logical_and(tip4 >= low_thresh, tip4 <= up_thresh)
        wheretipn =np.logical_and(notip >= low_thresh, notip <= up_thresh)
        
        ntip1 = np.where(wheretip1 == True)
        tip1_area = np.sum(tot_area1[ntip1])
        ntip3 = np.where(wheretip3 == True)
        tip3_area = np.sum(tot_area3[ntip3])
        ntip4 = np.where(wheretip4 == True)
        tip4_area = np.sum(tot_area4[ntip4])
        ntipn = np.where(wheretipn == True)
        ntip_area = np.sum(tot_arean[ntipn])
        
        if np.size(tipped) + np.size(ntip1) + np.size(ntip3)+ np.size(ntip3) + np.size(ntipn) ==0:
            tip_count[i] = np.nan
        else:
            tip_count[i] = tipped_area/(tipped_area + tip1_area + tip3_area+ tip4_area + ntip_area)
    
    colors = ['tab:blue', 'tab:green', 'tab:purple', 'darkred', 'dimgray', 'goldenrod', 'navy']     
    ax[7].plot(thresholds-0.1, tip_count*100, color=colors[q])    
    
    
    print ('done round '+str(q))
    
fig.text(0.5, -0.01, 'Sensitivity of seasonal cycle amplitude to global warming ($KK^{-1}$)', ha='center')
fig.text(-0.01, 0.5, '% of NSA region that suffers an AS', va='center', rotation='vertical')
ax[7].set_ylabel('Risk of an AS (%)')
ax[7].set_xlim([-0.6, 2.0])
ax[7].text(0.88, 0.88, '('+alphabet[8]+')', transform=ax[8].transAxes)
ax[7].set_xticks(np.arange(-0.5, 2.5, 0.5))
ax[7].set_yticks(np.arange(0,125,25))
ax[7].legend(models2, loc='best', bbox_to_anchor=(2, -0.5), ncol = 4)
plt.xticks(np.arange(-1, 1.1, 0.2))
ax[7].text(0.88, 0.88, '('+alphabet[7]+')', transform=ax[7].transAxes)

plt.yticks(np.arange(0, 1.1, 0.1))
plt.show()
    
savepath = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/Figures/Paper/Supplementary/Bootstap/'
filename = savepath+'boostrap_fig.svg'
fig.savefig(filename, format = 'svg', dpi=300, bbox_inches='tight')
        




