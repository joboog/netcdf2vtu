# -*- coding: utf-8 -*-

from netCDF4 import Dataset

import netcdf2vtk.netcdf2vtk as n2v

class Mapper(object):
    """A class for convenient mapping of a NetCDF to VTK."""

    def __init__(
        self,
        nc_path,
        vtu_path,
        data_var_names,
        map_func_type,
        nc_crs,
        vtu_crs,
        nullvalue,
        **kwargs_nc
    ):
        self.nc = dict(path = nc_path,
                        crs = nc_crs,
                        data_names = data_var_names)
        self.in_vtu = dict(path = vtu_path,
                        crs = vtu_crs)

        self.map_func_type = map_func_type
        self.nullvalue = nullvalue

        # read source nc
        dset = Dataset(self.nc["path"], mode = "r", format = "NETCDF4",
                        **kwargs_nc)
        self.nc.update({"dataset" : dset})

        # read dst vtu
        dst_vtu = n2v.read_vtu(self.in_vtu["path"])
        self.in_vtu.update({"vtu" : dst_vtu})


    def get_nc_variables(self):

        self.nc["dataset"].variables


    def set_nc_coord_names(self,
                            lat_name,
                            lon_name,
                            time_name = None):

        self.nc.update({"coord_names" : {"lat_name" : lat_name,
                                        "lon_name" : lon_name,
                                        "time_name" : time_name}})


    def set_nc_data_names(self, var_names):

        self.nc["data_names"] = var_names


    def read_nc_coords(self):

        lat_name, lon_name, time_name = self.nc["coord_names"].values()

        lat_dat = self.nc["dataset"].variables[lat_name][:].filled()
        lon_dat = self.nc["dataset"].variables[lon_name][:].filled()

        if time_name is not None:
            time = self.nc["dataset"].variables[time_name][:].filled()
        else:
            time = None

        coords = {"lat": lat_dat, "lon": lon_dat, "time" : time}
        self.nc.update({"coords" : coords})


    def read_nc_data(self):

        nc_data = n2v.get_nc_data(self.nc["dataset"],
                                  self.nc["data_names"])
        self.nc.update({"data" : nc_data})


    def interpolate(self):

        self.vtp = n2v.create_vtp(
                    self.nc["coords"]["lon"],
                    self.nc["coords"]["lat"],
                    self.nc["crs"],
                    self.in_vtu["crs"])

        n2v.nc_data_to_vtp(self.vtp,
                           self.nc["data"],
                           self.nc["coords"]["time"])

        print("NetCDF converted to VTP.")


        self.out_vtu = n2v.interpolate_vtp_data_on_vtu(
                        self.vtp,
                        self.in_vtu["vtu"],
                        self.map_func_type,
                        self.nullvalue)

        print("Data from VTP interpolated to VTU.")


    def write_out_vtu(self, path):
        # output updated ogs-mesh
        n2v.write_vtu(self.out_vtu, path)
        print("New VTU mesh written to disk.")


    def write_vtp(self, path):

        n2v.write_vtp(self.vtp, path)


    def map(self,
            out_path,
            lat_name,
            lon_name,
            time_name = None):
        """Map data from NetCDF to VTU."""

        self.set_nc_coord_names(lat_name,
                                lon_name,
                                time_name)
        self.read_nc_coords()
        self.read_nc_data()
        self.interpolate()
        self.write_out_vtu(out_path)
