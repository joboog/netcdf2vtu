# -*- coding: utf-8 -*-

# %%
from netCDF4 import Dataset
import netcdf2vtk.netcdf2vtk as nv2
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa

# %%
# extract coordinate extension from Selke_3D_Top.vtu
selke_vtu = nv2.read_ogs_vtu("data/Selke_3D_Top.vtu")

# %% get data objects
selke = dsa.WrapDataObject(selke_vtu)

# get point array
selke_points = selke.Points
xmin = selke_points[:,0].min()
xmax = selke_points[:,0].max()
ymin = selke_points[:,1].min()
ymax = selke_points[:,1].max()

#%% transfrom to crs of nc
from pyproj import Transformer

nc_crs = "EPSG:4326"
selke_crs = "EPSG:5684"

transf = Transformer.from_crs(selke_crs, nc_crs,                        always_xy=True)

lon_min, lat_min = transf.transform(xmin, ymin)
lon_max, lat_max = transf.transform(xmax, ymax)

#%% get nc
src_nc = Dataset("data/mhm.nc", "r", format="NETCDF4")

# %% get coordinate variables nc
nc_lons = src_nc.variables["lon"][:]
nc_lats = src_nc.variables["lat"][:]

# diff in lat
nc_lats[0,0]-nc_lats[1,0]
# %% diff in lon
nc_lons[0,1]-nc_lons[0,0]

# %% compute extension for cropping
nc_diff = 0.2
nc_lon_max = lon_max + nc_diff
nc_lon_min = lon_min - nc_diff
nc_lat_max = lat_max + nc_diff
nc_lat_min = lat_min - nc_diff

# %% get col index of nc_lon_min in nc_lons
nc_lon_min_bol = nc_lons < nc_lon_min
true_cols_lon = [not nc_lon_min_bol[:,i] for i in range(nc_lon_min_bol.shape[1])]
lon_min_idcol = min([i for i,val in enumerate(true_cols_lon) if val])

# %% get col index of nc_lon_min in nc_lons
nc_lon_min_bol = nc_lons < nc_lon_min
true_cols_lon = [all(nc_lon_min_bol[:,i]) for i in range(nc_lon_min_bol.shape[1])]
lon_min_idcol = max([i for i,val in enumerate(true_cols_lon) if val])

# %% get col index of nc_lon_max in nc_lons
nc_lon_bol = nc_lons > nc_lon_max
true_cols_lon = [all(nc_lon_bol[:,i]) for i in range(nc_lon_bol.shape[1])]
lon_max_idcol = min([i for i,val in enumerate(true_cols_lon) if val])

# %% get row index of nc_lat_min in nc_lats
nc_lat_bol = nc_lats < nc_lat_min
true_cols_lat = [all(nc_lat_bol[i,:]) for i in range(nc_lat_bol.shape[0])]
lat_min_idcol = min([i for i,val in enumerate(true_cols_lat) if val])

# %% get row index of nc_lat_max in nc_lats
nc_lat_bol = nc_lats > nc_lat_max
true_cols_lat = [all(nc_lat_bol[i,:]) for i in range(nc_lat_bol.shape[0])]
lat_max_idcol = max([i for i,val in enumerate(true_cols_lat) if val])

# %% subset lons and lats of src_nc
nc_lons_sub = nc_lons[lat_max_idcol:lat_min_idcol,
                    lon_min_idcol:lon_max_idcol]

nc_lats_sub = nc_lats[lat_max_idcol:lat_min_idcol,
                    lon_min_idcol:lon_max_idcol]
# %% get swc variables
nc_swcl01 = src_nc.variables["SWC_L01"][0:1]
nc_swcl02 = src_nc.variables["SWC_L02"][0:1]


nc_swcl01_sub = nc_swcl01[:,lat_max_idcol:lat_min_idcol,
                    lon_min_idcol:lon_max_idcol]
nc_swcl02_sub = nc_swcl02[:,lat_max_idcol:lat_min_idcol,
                    lon_min_idcol:lon_max_idcol]
# %%



# %%
# new create new nc file with subsets of data/mhm.nc
rootgrp = Dataset("mhm_test.nc", "w", format="NETCDF4")

# %%
# create dimensions
x = rootgrp.createDimension("x", nc_lons_sub.shape[1])
y = rootgrp.createDimension("y", nc_lons_sub.shape[0])
time = rootgrp.createDimension("time", 2)

# %%
# create variables
lon = rootgrp.createVariable("lon","f4", ("y","x"))
lat = rootgrp.createVariable("lat","f4", ("y","x"))
time = rootgrp.createVariable("time","f8", ("time"))
SWC_L01 = rootgrp.createVariable("SWC_L01","f8", ("time", "y", "x"), fill_value = -9999.0)
SWC_L02 = rootgrp.createVariable("SWC_L02","f8", ("time", "y", "x"), fill_value = -9999.0)


#%%
# add attributes
lon.setncattr_string("standard_name", "longitude")
lon.setncattr_string("long_name", "longitude")
lon.setncattr_string("units", "degrees_east")
lon.setncattr_string("_CoordinateAxisType", "Lon")


#%%
# add attributes
lat.setncattr_string("standard_name", "latitude")
lat.setncattr_string("long_name", "latitude")
lat.setncattr_string("units", "degrees_east")
lat.setncattr_string("_CoordinateAxisType", "Lat")

#%%
# add attributes
time.setncattr_string("standard_name", "time")
time.setncattr_string("long_name", "time")
time.setncattr_string("units", "hours since 2016-1-1 00:00:00")
time.setncattr_string("calendar", "standard")
time.setncattr_string("axis", "T")

#%%
# add values
lon[:] = nc_lons_sub
lat[:] = nc_lats_sub
time[:] = src_nc.variables["time"][0:2]

# %%
# add attributes
SWC_L01.setncattr_string("long_name", "soil water content of soil layer    1")
SWC_L01.setncattr_string("coordinates", "lat lon")
SWC_L01.setncattr("missing_value", -9999.)
SWC_L01.setncattr_string("unit", "mm")

# %%
# add attributes
SWC_L02.setncattr_string("long_name", "soil water content of soil layer    2")
SWC_L02.setncattr_string("coordinates", "lat lon")
SWC_L02.setncattr("missing_value", -9999.)
SWC_L02.setncattr_string("unit", "mm")

# %%
# add data to variables
SWC_L01[:] = nc_swcl01_sub
SWC_L02[:] = nc_swcl02_sub

# %%
rootgrp.close()



#%% now do the same with data/precip.nc
src_nc = Dataset("data/precip.nc", "r", format="NETCDF4")

# %% get coordinate variables nc
nc_lons = src_nc.variables["nlon"][:]
nc_lats = src_nc.variables["nlat"][:]

# unfortunately the extend of precip.nc does not cover the selke catchment
# lets create a similar nc file from dummy data


# %%
# new create new nc file with subsets of data/precip.nc
rootgrp = Dataset("precip_test.nc", "w", format="NETCDF4")

#%% create now coordiantes
nc_lons = np.arange(11,12,0.25)
nc_lats = np.arange(51.5,52.5,0.25)

# %%
# create dimensions
nlon = rootgrp.createDimension("nlon", len(nc_lons))
nlat = rootgrp.createDimension("nlat", len(nc_lats))

# %%
# create variables
nlon = rootgrp.createVariable("nlon","f4", ("nlon"))
nlat = rootgrp.createVariable("nlat","f4", ("nlat"))

precip = rootgrp.createVariable("precipitation","f8", ("nlon", "nlat"), fill_value = -9999.0)
relerr = rootgrp.createVariable("relativeError","f8", ("nlon", "nlat"), fill_value = -9999.0)


#%%
# add attributes
nlon.setncattr_string("long_name", "longitude")
nlon.setncattr_string("standard_name", "longitude")
nlon.setncattr_string("units", "degrees_east")

#%%
# add attributes
nlat.setncattr_string("standard_name", "latitude")
nlat.setncattr_string("long_name", "latitude")
nlat.setncattr_string("units", "degrees_north")

#%%
# add attributes
precip.setncattr_string("units", "mm/hr")
precip.setncattr_string("coordinates", "nlon nlat")
relerr.setncattr_string("units", "mm/hr")
relerr.setncattr_string("coordinates", "nlon nlat")


#%% add data
nlon[:] = nc_lons
nlat[:] = nc_lats

precip[:] = np.ma.MaskedArray(
                data = (np.random.random_sample(size=16)*0.05).reshape(4,4),
                mask = np.repeat(False, 16).reshape(4,4),
                fill_value = -9999.0)

relerr[:] = np.ma.MaskedArray(
                data = (np.random.random_sample(size=16)*0.0025).reshape(4,4),
                mask = np.repeat(False, 16).reshape(4,4),
                fill_value = -9999.0)

#%% done
rootgrp.close()













# %%
