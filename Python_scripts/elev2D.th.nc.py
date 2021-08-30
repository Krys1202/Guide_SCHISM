# @author: Krys1202
# Script to create elev2D.th.nc for SCHISM
# input: hgrid.ll
# output: elev2D.th.nc

from netCDF4 import Dataset,num2date
from datetime import datetime
import numpy as np
import xarray as xr

# get lon/lat of open boundaries nodes from hgrid.ll
f=open('hgrid.ll','r')
data=f.readlines()

num_nodes_and_elements = data[1]
num_nodes = num_nodes_and_elements.split(" ")[1]
num_elements= num_nodes_and_elements.split(" ")[0]

nodes_start = 2
node_lines = data[nodes_start:nodes_start+int(num_nodes)]

elements_start = nodes_start + int(num_nodes) 
element_lines = data[elements_start:elements_start + int(num_elements)]

metadata1_start = elements_start + int(num_elements) 
metadata1_line = data[metadata1_start:metadata1_start + 3]
metadata_obnb=data[metadata1_start+1]
numnodes_obnb=metadata_obnb.split(" ")[0] #to have total number of nodes in open bnd
metadata_obnb1=data[metadata1_start+2]
numnodes_obnb1=metadata_obnb1.split(" ")[0] 

obnodes1_start =metadata1_start+ 3
obnnode1_lines=data[obnodes1_start:obnodes1_start+int(numnodes_obnb1)]

metadata2_start = obnodes1_start + int(numnodes_obnb1)
metadata_obnb2=data[metadata2_start]
numnodes_obnb2=metadata_obnb2.split(" ")[0]

obnodes2_start =metadata2_start+ 1
obnnode2_lines=data[obnodes2_start:obnodes2_start+int(numnodes_obnb2)]
openbnd_lines=obnnode1_lines+obnnode2_lines #all the lines for all nodes in open bnd

#get the coordinates corresponding to the open bnd nodes
bnd_coord=[]
nodeLon=[]
nodeLat=[]
nodeLont=[]
nodeLatt=[]
for i in range(len(openbnd_lines)):
        node_index = int(openbnd_lines[i])
        node= node_lines[node_index+2]
        #node_num=node[0:5]
        node_lont=node.split(" ")[5]
        node_latt=node.split(" ")[10]
        #append lon, lat
        nodeLont.append(node_lont)
        nodeLatt.append(node_latt)

bnd_coord= np.column_stack([nodeLont, nodeLatt])

#nOpenBnbNodes=int(numnodes_obnb)

##create new netcdf elev2D.th.nc
filename = 'elev2Dmarch.th.nc'
nc_elev = Dataset(filename, 'w', format='NETCDF4')
#nc_elev = Dataset(filename, 'w', format='NETCDF3_CLASSIC')
nc_elev.setncatts({"Conventions": "CF-1.0"})
nc_elev.createDimension('time', None)
nc_elev.createDimension('nComponents', 1) 
nc_elev.createDimension('nLevels', 1) 
nc_elev.createDimension('nOpenBndNodes', len(bnd_coord)) #nOpenBnbNodes
nc_elev.createDimension('one', 1) 

one=1
time_var = nc_elev.createVariable('time', np.float32, ('time',))
time_var.units = 's'
time_var.long_name = 'simulation time in seconds'

timestep_var = nc_elev.createVariable('time_step', np.float32, 'one')
timestep_var.units = 's'
timestep_var.long_name = 'time step in seconds'
timestep_var=86400

elev_var = nc_elev.createVariable('time_series',np.float32,('time','nOpenBndNodes',
                             'nLevels', 'nComponents'),)
elev_var.long_name = 'elev_nontidal'
elev_var.units='m'

#Fill the new nc file with Mercator data ( ex 1month  01/07/2018)
d  = Dataset('merc-all-201803-1s.nc') #modified from original version with 
#cdo to have time in seconds
tvar = 'time'
t = d[tvar][:]
t_unit = d.variables[tvar].units
#print(t_unit)
tvals = num2date(t,units = t_unit)
str_t = [i.strftime("%Y%m%d %H%M") for i in tvals] # to display dates as string
datetime_t = [datetime.strptime(i,"%Y%m%d %H%M") for i in str_t]
time_var[:]= t


# need to find the values for lon/lat of openbnb =bnd_coord
nc=xr.open_dataset('merc-all-201803-1s.nc') 
subset=nc.sel(longitude=bnd_coord[:,0], latitude=bnd_coord[:,1], method='nearest')
elev_test=subset.zos.values

elev_test2=np.zeros([len(t),len(bnd_coord)])

for j in range(len(t)):
    for i in range(len(bnd_coord)):
      elev_test2[j,i] = elev_test[j,i,i]

#put nan values=0
elev_nonan=np.nan_to_num(elev_test2)
elev_var[:,:,:,:] = elev_nonan

d.close()
nc_elev.close()
