#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from netcdf2vtk.netcdf2vtk import Mapper

# def variables ------------------------------------------------
nc_path = "data/ex2_4.nc" # path of netcdf file
vtu_path = "data/ogs.vtu" # path of ogs-mesh
vtu_new_path = "ex4_new.vtu" # path of updated ogs-mesh
data_var_names = ["precipitation", "relativeError"] # names of the netcdf data variables
map_func_type = 1 # def mapping func type 1: Voronoi, 2:Gaussian, 3:Shepard
nc_crs = "EPSG:4326"   # coordinate system of netcdf file
vtu_crs = "EPSG:5684"   # coordinate sysetm of ogs-mesh file
nullvalue = -9999

# run script ------------------------------------------
# create mapper object
mapper = Mapper(nc_path,
                    vtu_path,
                    data_var_names,
                    map_func_type,
                    nc_crs,
                    vtu_crs,
                    nullvalue)

# set the names of coordinates variables in nc file
mapper.set_nc_coord_names(lat_name = "nlat",
                        lon_name = "nlon")

# extract coordinate and variable data from nc file
mapper.read_nc_coords()
mapper.read_nc_data()

# interpolate data from nc to vtu
mapper.interpolate()

# output the VTU with the mapped data
mapper.write_out_vtu(path = vtu_new_path)


# for convenience, all the above can be within two statements
mapper2 = Mapper(nc_path,
                    vtu_path,
                    data_var_names,
                    map_func_type,
                    nc_crs,
                    vtu_crs,
                    nullvalue)

mapper2.map(out_path = "ex4_mapper2_new_vtu.vtu",
            lat_name = "nlat",
            lon_name = "nlon")