
# coding: utf-8

# In[1]:

#get_ipython().magic(u'matplotlib inline')

import numpy as np
import pandas as pd
from math import *
import matplotlib.pyplot as plt
import scipy.io as sio

from mpl_toolkits.basemap import Basemap, cm


# In[2]:

# load in the NPP data

fname = '/Users/sclayton/Documents/MIT_work/ostreo/satellite_data/vgpm.2009289.all.xyz'
result = pd.read_csv(fname, sep = ' ')
print(result.columns)


# In[3]:

lon = (result['lon']).values # longitude
lat = (result['lat']).values # latitude
npp = (result['value']).values # npp

# reshape the data to a 2160x4320 grid
ny = 2160
nx = 4320

LON = lon.reshape(ny,nx)
LAT = lat.reshape(ny,nx)
NPP = npp.reshape(ny,nx)



# In[ ]:

# import location data for cruise track and stations
cfname = '/Volumes/sclayton/kuroshio/AVHRR_9km/xbt_ctd_coords.mat'
coords = sio.loadmat(cfname)

botlat = coords['bot_lat']
botlon = coords['bot_lon']
xbtlon = coords['xbt_lon']
xbtlat = coords['xbt_lat']

print xbtlat


# In[ ]:

# make a plot of the data

fig1 = plt.figure(1, figsize=(20,12))

m = Basemap(width=1500000,height=1200000,projection='lcc',
            resolution='c',lat_1=30.,lat_2=50,lat_0=37,lon_0=145.)

m.drawcoastlines()
m.drawmapboundary(fill_color='white')
m.fillcontinents(color='grey',lake_color='white')

im1 = m.pcolor(LON,LAT,NPP, edgecolor = 'none' , cmap=plt.cm.jet, latlon=True, vmin= 0, vmax = 2000)
m.plot(xbtlon,xbtlat)

m.drawparallels(np.arange(30,50,5),labels=[1,1,0,1])
m.drawmeridians(np.arange(138,152,2),labels=[1,1,0,1])

cb = plt.colorbar(im1, shrink=.92)
cb.set_label(label='NPP (mg/C/m^2/day)',fontsize=16)

plt.tick_params(axis='both', which='major', labelsize=16)

plt.show()


# In[ ]:



