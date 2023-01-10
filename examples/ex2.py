#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: joboog
"""
import netcdf2vtk.netcdf2vtk as n2v
from netCDF4 import Dataset

# def variables ------------------------------------------------
src_nc_path = "data/precip.nc" # path of netcdf file
ogs_vtu_path = "data/Selke_3D_Top.vtu" # path of ogs-mesh
ogs_vtu_new_path = "ex2_new2.vtu" # path of updated ogs-mesh
data_var_names = ["precipitation", "relativeError"] # names of the netcdf data variables
map_func_type = 1 # def mapping func type 1: Voronoi, 2:Gaussian, 3:Shepard
src_nc_crs = "EPSG:4326"   # coordinate system of netcdf file
ogs_vtu_crs = "EPSG:5684"   # coordinate sysetm of ogs-mesh file
nullvalue = -9999

# run script -----------------------------------------------------
# import netcdf
src_nc = Dataset(src_nc_path, mode = "r", format = "NETCDF4")

# get spatial and temporal coordinates
lat_dat = src_nc.variables["nlat"][:].filled()
lon_dat = src_nc.variables["nlon"][:].filled()
#time = src_nc.variables["time"][:].filled()

# extract netcdf data and transform to vtkPolyData
src_vars = n2v.get_src_nc_data(src_nc, data_var_names)
src_poly = n2v.init_src_poly(lon_dat, lat_dat, src_nc_crs, ogs_vtu_crs)
n2v.add_nc_data_to_src_poly(src_poly, src_vars)
print("finish")

# import ogs-mesh and map netcdf data on
ogs_vtu = n2v.read_ogs_vtu(ogs_vtu_path)
print("ogs-mesh read")
ogs_vtu_new = n2v.map_data_on_ogs_vtu(src_poly, ogs_vtu, map_func_type, nullvalue)
print("data mapped")
# output updated ogs-mesh
n2v.write_mapped_ogs_vtu(ogs_vtu_new, ogs_vtu_new_path)
print("ogs-mesh updated")