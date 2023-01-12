#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import netcdf2vtk.mapper as nvmp

# def variables ------------------------------------------------
src_nc_path = "data/ex2_4.nc" # path of netcdf file
ogs_vtu_path = "data/ogs.vtu" # path of ogs-mesh
ogs_vtu_new_path = "ex4_new.vtu" # path of updated ogs-mesh
data_var_names = ["precipitation", "relativeError"] # names of the netcdf data variables
map_func_type = 1 # def mapping func type 1: Voronoi, 2:Gaussian, 3:Shepard
src_nc_crs = "EPSG:4326"   # coordinate system of netcdf file
ogs_vtu_crs = "EPSG:5684"   # coordinate sysetm of ogs-mesh file
nullvalue = -9999

# run script ------------------------------------------
# create mapper object
mapper = nvmp.Mapper(src_nc_path,
                    ogs_vtu_path,
                    data_var_names,
                    map_func_type,
                    src_nc_crs,
                    ogs_vtu_crs,
                    nullvalue)

# extract data from src file
mapper.get_scr_nc_data()

# initialize a VTP object from coordinates in the netcdf file
mapper.get_src_coords(lat_name = "nlat", lon_name = "nlon")
mapper.init_src_poly()

# map netcdf data via the created VTP on the destination vtu
mapper.map_data()

# output the VTU with the mapped data
mapper.write_out_vtu(path = ogs_vtu_new_path)


# for convenience, all the above can be within two statements
mapper2 = nvmp.Mapper(src_nc_path,
                    ogs_vtu_path,
                    data_var_names,
                    map_func_type,
                    src_nc_crs,
                    ogs_vtu_crs,
                    nullvalue)

mapper2.map(out_path="ex4_mapper2_new_vtu.vtu",
            lat_name="nlat",
            lon_name="nlon")