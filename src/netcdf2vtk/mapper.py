# -*- coding: utf-8 -*-

from netCDF4 import Dataset

import netcdf2vtk.netcdf2vtk as n2v

class Mapper(object):
    """A class for convenient mapping of a NetCDF to VTK."""

    def __init__(
        self,
        src_nc_path,
        dst_vtu_path,
        data_var_names,
        map_func_type,
        src_crs,
        dst_crs,
        nullvalue
    ):
        self.src_nc_path = src_nc_path
        self.dst_vtu_path = dst_vtu_path
        self.data_var_names = data_var_names
        self.map_func_type = map_func_type
        self.src_crs = src_crs
        self.dst_crs = dst_crs
        self.nullvalue = nullvalue

        # read source nc
        self.src_nc = Dataset(self.src_nc_path, mode = "r", format = "NETCDF4")
        # read dst vtu
        self.dst_vtu = n2v.read_ogs_vtu(self.dst_vtu_path)


    def get_src_coords(self,
                        lat_name,
                        lon_name,
                        time_name = None):

        lat_dat = self.src_nc.variables[lat_name][:].filled()
        lon_dat = self.src_nc.variables[lon_name][:].filled()

        if time_name is not None:
            time = self.src_nc.variables[time_name][:].filled()
        else:
            time = None

        self.src_coords = {"lat": lat_dat, "lon": lon_dat, "time" : time}


    def set_data_var_names(self, var_names):

        self.data_var_names = var_names


    def get_scr_nc_data(self):

        self.src_data = n2v.get_src_nc_data(self.src_nc, self.data_var_names)


    def init_src_poly(self):

        self.src_poly = n2v.init_src_poly(
                            self.src_coords["lon"],
                            self.src_coords["lat"],
                            self.src_crs,
                            self.dst_crs)

        n2v.add_nc_data_to_src_poly(self.src_poly,
                                    self.src_data,
                                    self.src_coords["time"])

        print("VTP created from NetCDF.")


    def map_data(self):

        self.out_vtu = n2v.map_data_on_ogs_vtu(
                        self.src_poly, self.dst_vtu, self.map_func_type,
                        self.nullvalue)
        print("Data mapped.")


    def write_out_vtu(self, path):
        # output updated ogs-mesh
        n2v.write_mapped_ogs_vtu(self.out_vtu, path)
        print("New VTU mesh written to disk.")


    def write_src_vtp(self, path):
        n2v.write_src_poly(self.src_poly, path)


    def map(self,
            out_path,
            lat_name,
            lon_name,
            time_name = None):
        """Map data from NetCDF to VTU."""
        self.get_src_coords(lat_name, lon_name, time_name)
        self.get_scr_nc_data()
        self.init_src_poly()
        self.map_data()
        self.write_out_vtu(out_path)
