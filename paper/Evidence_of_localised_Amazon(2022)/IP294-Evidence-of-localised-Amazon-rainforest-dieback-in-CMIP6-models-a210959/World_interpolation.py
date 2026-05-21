#- coding: utf-8 -*-
"""
Created on Mon Dec 16 09:01:27 2019

@author: pdlr201

Script to interpolate annual 1pctCO2 CMIP6 data onto 1 degree world grid
"""

import iris
import numpy as np
import numpy.ma as ma
from iris.cube import Cube


def shiftlon(lon,lon_0):
    """returns original sequence of longitudes (in degrees) recentered
    in the interval [lon_0-180,lon_0+180]"""
    lon_shift = np.asarray(lon)
    lon_shift = np.where(lon_shift > lon_0+180, lon_shift-360 ,lon_shift)
    lon_shift = np.where(lon_shift < lon_0-180, lon_shift+360 ,lon_shift)
    itemindex = len(lon)-np.where(lon_shift[0:-1]-lon_shift[1:] >= 180)[0]
    return (np.roll(lon_shift,itemindex-1), itemindex-1)

# You shouldn't need to change the below entries at least for the 3 models you
# are currently using
temp_res = 'Lmon'
experiment = '1pctCO2'
variant_id = 'r1i1p1f1'

################# Specify variable of interest ######################
var = 'cVeg' # or treeFrac

################# Specify model here ################################## 
model = 'GFDL-ESM4'


# Specific entries for the 3 models initially used
if model == 'EC-Earth3-Veg':

    grid_label = '_gr'
    date_ranges = []
    for k in np.arange(151):
        date_ranges.extend([str(1850+k)+'01-'+str(1850+k)+'12'])
    total_years = 151
    date_range2 = '1850-2000'

elif model == 'GFDL-ESM4':
    
    grid_label = '_gr1'
    date_ranges = ['000101-010012', '010101-015012']
    total_years = 150
    date_range2 = '1850-1999'
    
elif model == 'MPI-ESM1-2-LR':
    
    grid_label = '_gn'
    date_ranges = []
    for k in np.arange(8):
        date_ranges.extend([str(1850+k*20)+'01-'+str(1850+((k+1)*20)-1)+'12'])
    date_ranges.extend(['201001-201412'])
    total_years = 165
    date_range2 = '1850-2014'

elif model == 'NorCPM1':
    grid_label = '_gn'
    date_ranges = ['000102-016412']
    total_years = 164
    date_range2 = '1850-2013'

elif model == 'TaiESM1':
    grid_label = '_gn'
    date_ranges = ['000102-015012']
    total_years = 150
    date_range2 = '1850-2000'
    
elif model == 'SAM0-UNICON':
    grid_label = '_gn'
    date_ranges = []
    for k in np.arange(15):
        date_ranges.extend([str(1850+(k*10))+'01-'+str(1859+(k*10))+'12'])
    total_years = 150
    date_range2 = '1850-1999'
    
elif model == 'UKESM1-0-LL':
    grid_label = '_gn'
    date_ranges = ['185001-194912', '195001-199912']
    total_years = 150
    date_range2 = '1850-2000'
    
######################## Specify path to data ##########################
# (I have 2 directories a raw data directory and a processed data directory)
path = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/'+var+'/'+experiment

######################## Specify directory of raw data ##################
raw_data = '/Original_data/'

######################## Specify directory of processed data ##################
processed_data = '/Processed_data/'

############ Specify directory within processed data directory for the region interpolated #################
region = 'World'

# Define latitude and longitude (uses 1 degree by 1 degree grid) 
my_lats=['latitude',np.arange(-90,91,step=1)]
my_lons=['longitude',np.arange(-180,181,step=1)]

# Months in the year
freq = 12

# Initialise counters
count = 0
sub_total_years = 0

# Loop through data files
for i in range(len(date_ranges)):
    
    # File name of data
    fname = path+raw_data+var+'_'+temp_res+'_'+model+'_'+experiment+'_'+variant_id+grid_label+'_'+date_ranges[i]+'.nc'

    # Load file contents into an iris cube
    x = iris.load_cube(fname)
    
    # Extract data values
    data = x.data

    # Determine lengths of data dimensions
    nt = int(len(data[:,0,0]))
    ny = int(len(data[0,:,0]))
    nx = int(len(data[0,0,:]))
    
    years = int((nt+1)/freq)
    
    # On first loop through create a new empty cube with original coordinates
    if count == 0:
    
        new_coord = iris.coords.DimCoord(range(total_years), long_name='time', units=1)
        coord1 = x.coord('latitude')
        coord2 = x.coord('longitude')
        
        cube = Cube(ma.zeros((total_years,ny,nx),np.float32),dim_coords_and_dims=[(new_coord,0),(coord1,1),(coord2,2)])
    
    count2 = 0
    
    # Stack all data into 1 cube 
    for j in range(sub_total_years,sub_total_years+years):

        cube.data[j,:,:] = np.mean(data[count2*freq:(count2+1)*freq,:,:], axis=0)
        
        count2 = count2 + 1
  
    sub_total_years = sub_total_years + years
    count = count+1

# Centre coordinates about longitude 0
cube.coord('longitude').points, rollindex = shiftlon(cube.coord('longitude').points, 0)

# Roll data to be consistent with longitude 0
cube.data = np.roll(cube.data, rollindex, axis=2)

# Interpolate data onto new coordinate system
region_cube = cube.interpolate([my_lats, my_lons],iris.analysis.Linear(extrapolation_mode='nan'))

# Ensure new cube has correct name
region_cube.rename(var)

# New file name for data
outname = path+processed_data+region+'/'+var+'_'+model+'_'+experiment+'_'+variant_id+'_'+date_range2+'_'+region+'.nc'

# Save data to file name
iris.save(region_cube, outname, unlimited_dimensions=[])
