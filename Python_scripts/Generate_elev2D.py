### THis script is for generating non-tidal forcing from both Mercator and HYCOM datasets
### The following is an example of generating elev2D.th.nc from Mercator forcing. WE can use the same script for HYCOM as well by just changing data variable names

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
from pyschism.mesh.vgrid import Vgrid
from pyschism.mesh.base import Nodes, Elements
from pyschism.mesh.hgrid import Hgrid
import xarray as xr

#### EXTRACTING OPEN BOUNDARY NODES AND ITS CORDINATES ################################################
# some of these lines are adapted from PYSCHISM
hgrid = Hgrid.open('./hgrid.gr3', crs='epsg:4326')
#vgrid = Vgrid.open('./vgrid.in')
gdf=hgrid.boundaries.open.copy()
opbd=[]
for ibnd in [0,1]: 
    opbd.extend(list(gdf.iloc[ibnd].indexes)) #extend adds a list or tuple at the end of a list while append adds element
    # stach one openboundary after other

blon = hgrid.coords[opbd,0] # all longitudes for open boundary nodes
blat = hgrid.coords[opbd,1] # all latitudes for open boundary nodes

## open the Meractor file
Mercator_file=xr.open_dataset('./Mercator_GLORYS/cmems_mod_glo_phy_my_0.083deg_P1D-m_2015_second.nc') #daily forcing


#Mercator_file=xr.open_dataset('./cmems_mod_glo_phy_anfc_0.083deg_PT1H-m_1711515859809.nc') #hourly forcing
# Always check and set the start datetime same as the start datetime in param.nml otherwise the model will pick up wrong timestep foricng 
# because the timestep in the non-tidal forcing files (which we are genearting) has no timestamp in datetime. Timestep is only given as an index
# e.g like 0, 10800, ......
run_start_date='2015-11-15' # also provide timestep of model start if available in the forcing
run_end_date='2015-11-23' # should cover equal or more than model run length
Mercator_file=Mercator_file.sel(time=slice(run_start_date,run_end_date))



# Lets extract elevation values for our open boundary cordinates 'blon' and 'blat'
extracted_elevation_values = []
for t in Mercator_file.time:
    elevation_values_timestep = []
    for lat, lon in zip(blat,blon):
        elevation = Mercator_file.zos.sel(latitude=lat, longitude=lon, method='nearest', time=t).values #here I am using a simple nearest neighbor method, can also use other type
        elevation_values_timestep.append(elevation)
    extracted_elevation_values.append(elevation_values_timestep)

extracted_elevation_values=np.array(extracted_elevation_values)

### CHECK BY PLOTING THE EXTRACTED SEA SURFACE HEIGHT ###
#plt.figure(figsize=(6,4))
#plt.imshow(extracted_elevation_values, cmap='viridis', interpolation='nearest')
#plt.colorbar()  # Add colorbar
#plt.gca().set_aspect('auto') 
#plt.gca().invert_yaxis()
#plt.show()

# PLEASE NOTE THAT THERE WILL BE A DISTINCT DIFFERENCE IN COLOR (VALUES) ON THE PLOT BETWEEN THE TWO OPEN BOUNDARIES



# Ok now let's create a dummny xarray dataset and we will replace with the extracted values
# Define dimensions and data
n_timesteps = len(Mercator_file.time)
#time_interval_hours = 1 # hourly forcing, change this interval if other timestep is used
time_interval_hours = 24 #for daily, change this interval if other timestep is used
time_interval_seconds = time_interval_hours * 3600 # timesteps in second

# Assuming `bnd_coord` contains the boundary coordinates
nOpenBndNodes = len(opbd)
nLevels = 1
nComponents = 1
one=1
# Create time array
time = np.arange(n_timesteps) * time_interval_seconds
# Create xarray dataset
ds = xr.Dataset()

# Add dimensions
ds['time'] = time
ds['nOpenBndNodes'] = np.arange(nOpenBndNodes)
ds['nLevels'] = np.arange(nLevels)
ds['nComponents'] = np.arange(nComponents)
ds['one'] = np.arange(1)  # Adding 'one' dimension with value 1

# Add variables
ds['time_series'] = (('time', 'nOpenBndNodes', 'nLevels', 'nComponents'), np.random.rand(n_timesteps, nOpenBndNodes, nLevels, nComponents))
ds['time'].attrs['units'] = 's'
ds['time'].attrs['long_name'] = 'simulation time in seconds'

ds['time_step'] = xr.DataArray(np.array([time_interval_seconds]), dims='one')
ds['time_step'].attrs['units'] = 's'
ds['time_step'].attrs['long_name'] = 'model time step in seconds'

ds['time_series'].attrs['long_name'] = 'sea surface height' #GIVE ANY NAME YOU LIKE
ds['time_series'].attrs['units'] = 'm'

#### PADDING THE EXTRACTED VALUES IN OUR DUMMY DATASET 
##  CREATING A COPY is a choice, can skip as well.
ds_new = ds.copy()
for i, t in enumerate(ds['time']):
    ds_new['time_series'][i, :, :, :] = extracted_elevation_values[i].reshape((312, 1, 1)) #always check the shape otherwise, woring data may pad in

ds_new=ds_new.drop_vars(['one','nOpenBndNodes','nComponents','nLevels']) #not compulsory, 
#ds_new.time_series.plot() ## THis plot should exactly look like the plot of extracted elevation values, if not same, there is something wrong while assigning values

## Finally save it
## It is always better to remove FillValue attributes where models can't handle such values 
ds_new.to_netcdf('./elev2D.th_MERC_2015.nc',format='NETCDF4', unlimited_dims='time', encoding={'time_series': {'_FillValue':None,'dtype': 'float32'},
                                                                                             'time': {'_FillValue':None,'dtype': 'float32'},
                                                                                              'time_step': {'_FillValue':None,'dtype': 'float32'}})
