### THis script is for generating non-tidal forcing from both Mercator and HYCOM datasets
### The following is an example of generating uv3D.th.nc from Mercator forcing. WE can use the same script for HYCOM as well by just changing data variable names

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
from pyschism.mesh.vgrid import Vgrid
from pyschism.mesh.base import Nodes, Elements
from pyschism.mesh.hgrid import Hgrid
import xarray as xr

#### EXTRACTING OPEN BOUNDARY NODES AND ITS CORDINATES ################################################
# Some of these lines are adapted from PYSCHISM
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

print(f'{Mercator_file.time[0]}')

# Lets extract U AND V values for our open boundary cordinates 'blon' and 'blat'
extracted_u_values = []
for t in Mercator_file.time:
    u_values_timestep = []
    for lat, lon in zip(blat,blon):
        u = Mercator_file.uo.sel(latitude=lat, longitude=lon, method='nearest', time=t).values
        u_values_timestep.append(u)
    extracted_u_values.append(u_values_timestep)

extracted_u_values=np.array(extracted_u_values)

extracted_v_values = []
for t in Mercator_file.time:
    v_values_timestep = []
    for lat, lon in zip(blat,blon):
        v = Mercator_file.vo.sel(latitude=lat, longitude=lon, method='nearest', time=t).values
        v_values_timestep.append(v)
    extracted_v_values.append(v_values_timestep)

extracted_v_values=np.array(extracted_v_values)

### ALWAYS CHECK BY PLOTING THE EXTRACTED U AND V SPEED ###
#plt.figure(figsize=(6,4))
#plt.imshow(extracted_u_values, cmap='viridis', interpolation='nearest')
#plt.colorbar()  # Add colorbar
#plt.gca().set_aspect('auto') 
#plt.gca().invert_yaxis()
#plt.show()

# PLEASE NOTE THAT THERE WILL BE A DISTINCT DIFFERENCE IN COLOR (VALUES) ON THE PLOT BETWEEN THE TWO OPEN BOUNDARIES


# Ok now let's create a dummny xarray dataset and we will replace with the extracted values
### Similarly for u and v

# Define dimensions and data
n_timesteps = len(Mercator_file.time)
#time_interval_hours = 1 # hourly forcing, change this interval if other timestep is used
time_interval_hours = 24 # daily forcing, change this interval if other timestep is used
time_interval_seconds = time_interval_hours * 3600 # timesteps in second

# Assuming `bnd_coord` contains the boundary coordinates
nOpenBndNodes = len(opbd)
nLevels = 2
nComponents = 2
one=1
# Create time array
time = np.arange(n_timesteps) * time_interval_seconds
# Create xarray dataset
ds_uv = xr.Dataset()

# Add dimensions
ds_uv['time'] = time
ds_uv['nOpenBndNodes'] = np.arange(nOpenBndNodes)
ds_uv['nLevels'] = np.arange(nLevels)
ds_uv['nComponents'] = np.arange(nComponents)
ds_uv['one'] = np.arange(1)  # Adding 'one' dimension with value 1

# Add variables
ds_uv['time_series'] = (('time', 'nOpenBndNodes', 'nLevels', 'nComponents'), np.random.rand(n_timesteps, nOpenBndNodes, nLevels, nComponents))
ds_uv['time'].attrs['units'] = 's'
ds_uv['time'].attrs['long_name'] = 'simulation time in seconds'

ds_uv['time_step'] = xr.DataArray(np.array([time_interval_seconds]), dims='one')
ds_uv['time_step'].attrs['units'] = 's'
ds_uv['time_step'].attrs['long_name'] = 'model time step in seconds'

ds_uv['time_series'].attrs['long_name'] = 'non-tidal velocity'
ds_uv['time_series'].attrs['units'] = 'm/s'

#### PADDING THE EXTRACTED VALUES IN OUR DUMMY DATASET 
##  CREATING A COPY is a choice, can skip as well.
ds_uv_new = ds_uv.copy()
for i, t in enumerate(ds_uv['time']):
    ds_uv_new['time_series'][i,:,0,0] = extracted_u_values[i].squeeze()
    ds_uv_new['time_series'][i,:,1,1] = extracted_v_values[i].squeeze()

ds_uv_new=ds_uv_new.drop_vars(['one','nOpenBndNodes','nComponents','nLevels']) #not compulsory, 
#ds_uv_new.time_series.sel(nComponents=0,nLevels=0).plot(cmap='viridis') ## THis plot should exactly look like the plot of extracted elevation values, if not same, there is something wrong while assigning values

## Finally save it
## It is always better to remove FillValue attributes where models can't handle such values 
ds_uv_new.to_netcdf('./uv3D.th_MERC_2015.nc',format='NETCDF4', unlimited_dims='time', encoding={'time_series': {'_FillValue':None,'dtype': 'float32'},
                                                                                             'time': {'_FillValue':None,'dtype': 'float32'},
                                                                                              'time_step': {'_FillValue':None,'dtype': 'float32'}})
