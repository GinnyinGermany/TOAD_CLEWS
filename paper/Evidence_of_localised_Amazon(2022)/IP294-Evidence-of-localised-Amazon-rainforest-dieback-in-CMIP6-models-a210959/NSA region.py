# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 11:50:22 2021

@author: impy2
"""
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 10})

import numpy as np
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker
from cartopy.io import srtm

from cartopy.io import PostprocessedRasterSource, LocatedImage
from cartopy.io.srtm import SRTM3Source, SRTM1Source

def shade(located_elevations):
    """
    Given an array of elevations in a LocatedImage, add a relief (shadows) to
    give a realistic 3d appearance.

    """
    new_img = srtm.add_shading(located_elevations.image,
                               azimuth=135, altitude=15)
    return LocatedImage(new_img, located_elevations.extent)

lonstep = 10
latstep = 10

fig, ax = plt.subplots(figsize= (2.7, 2.0), subplot_kw=dict(projection=ccrs.PlateCarree()))
# ax = plt.axes(projection=ccrs.PlateCarree() )
ax.stock_img()

ax.add_feature(cfeature.COASTLINE, zorder=10)
NSA_verts = np.array([[-72,12], [-72,-8], [-50,-8], [-50,7.6], [-55,12], [-72,12]])
NWS_verts = np.array([[-75,12], [-83.4,2.2], [-83.4,-10], [-79,-15], [-72,-15], [-72,12], [-75,12]])
NES_verts = np.array([[-34,-20], [-50,-20], [-50,0], [-34,0], [-34, -20]])
SAM_verts = np.array([[-66.4,-20], [-72,-15], [-72,-8], [-50,-8], [-50,-20], [-66.4,-20]])
    
# NSA_lons = [-72, -72, -50, -50, -55, -72]
# NSA_lats = [12, -8, -8, 7.6, 12, 12]

# ax.set_ylim((-20, 15))
ax.text(-65, 0, 'NSA',
          horizontalalignment='left',
          transform=ccrs.Geodetic())
# ax.text(-78, -6, 'NWS',
#           horizontalalignment='left',
#           transform=ccrs.Geodetic(), rotation=90)
# ax.text(-46, -11, 'NES',
#           horizontalalignment='left',
#           transform=ccrs.Geodetic())
# ax.text(-65, -15, 'SAM',
#           horizontalalignment='left',
#           transform=ccrs.Geodetic())



lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()

ax.plot(NSA_verts[:,0], NSA_verts[:,1],
          color='black', linewidth=1,
          transform=ccrs.Geodetic())

# ax.plot(NWS_verts[:,0], NWS_verts[:,1],
#           color='black', linewidth=1,
#           transform=ccrs.Geodetic())

# ax.plot(NES_verts[:,0], NES_verts[:,1],
#           color='black', linewidth=1,
#           transform=ccrs.Geodetic())

# ax.plot(SAM_verts[:,0], SAM_verts[:,1],
#           color='black', linewidth=1,
#           transform=ccrs.Geodetic())


ax.gridlines(draw_labels=False, ylocs = np.arange(-20,20,10), linestyle='--')

ax.set_xticks(np.arange(-80,-30,20),crs=ccrs.PlateCarree())
ax.set_xticklabels(np.arange(-80,-30,20))
ax.set_yticks(np.arange(-20,20,10),crs=ccrs.PlateCarree())
ax.set_yticklabels(np.arange(-20,20,10))

lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()

ax.set_xticklabels(ax.get_xticks()) 

ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.grid(linewidth=0.5, color='black', linestyle='--') 
fig.text(0.85, 0.79, '(a)', ha='center')


ax.set_extent([-86,-33, -20, 13])


plt.show()

savepath = 'C:/Users/impy2/OneDrive/Documents/Uni Yr3/Tipping Points Project/Figures/Paper/Figure 2/'
filename = savepath+'NSA_region(2).svg'
fig.savefig(filename, format = 'svg', dpi=300, bbox_inches='tight') 

