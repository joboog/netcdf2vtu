# -*- coding: utf-8 -*-

from netCDF4 import Dataset
from numpy.testing import assert_array_almost_equal
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid, vtkPolyData
from os.path import isfile

import netcdf2vtk.netcdf2vtk as n2v


nc_file1 = "data/ex1_3.nc"
nc_file2 = "data/ex2_4.nc"

in_vtu_path = "data/ogs.vtu"
out_vtu_path = "temp.vtu"
out_vtp_path = "temp.vtp"
data_var_names = ["SWC_L01", "SWC_L02"]
map_func_type = 1
nc_crs = "EPSG:4326"
vtu_crs = "EPSG:5684"
nullvalue = -9999
coordnames = ["lat", "lon", "time"]

class TestMapper():

    def setup_method(self):
        self.nc_file = nc_file1
        self.nc_data_var = data_var_names
        self.vtu_in = in_vtu_path
        self.vtu_out = out_vtu_path
        self.vtp_out = out_vtp_path
        self.nc_crs = nc_crs
        self.vtu_crs = vtu_crs
        self.map_func_type = map_func_type
        self.nullval = nullvalue
        self.coordnames = coordnames

        self.dataset = Dataset(self.nc_file, mode = "r", format = "NETCDF4")
        self.vtu = n2v.read_vtu(self.vtu_in)

        self.mapper = n2v.Mapper(
            nc_path = self.nc_file,
            vtu_path = self.vtu_in,
            data_var_names = self.nc_data_var,
            map_func_type = self.map_func_type,
            nc_crs = self.nc_crs,
            vtu_crs = self.vtu_crs,
            nullvalue = self.nullval
        )

        self.lat_dat = self.dataset.variables[self.coordnames[0]][:].filled()
        self.lon_dat = self.dataset.variables[self.coordnames[1]][:].filled()
        self.time_dat = self.dataset.variables[self.coordnames[2]][:].filled()

        self.nc_data = n2v.get_nc_data(self.dataset, self.nc_data_var)


    def test_init_mapper(self):

        mapper = n2v.Mapper(
            nc_path = self.nc_file,
            vtu_path = self.vtu_in,
            data_var_names = self.nc_data_var,
            map_func_type = self.map_func_type,
            nc_crs = self.nc_crs,
            vtu_crs = self.vtu_crs,
            nullvalue = self.nullval
        )

        assert type(mapper) is n2v.Mapper

        assert mapper.nc["path"] == self.nc_file
        assert mapper.nc["crs"] == self.nc_crs
        assert mapper.nc["data_names"] == self.nc_data_var
        assert type(mapper.nc["dataset"]) == type(self.dataset)

        assert mapper.in_vtu["path"] == self.vtu_in
        assert mapper.in_vtu["crs"] == self.vtu_crs
        assert type(mapper.in_vtu["vtu"]) == type(self.vtu)

        assert mapper.map_func_type == self.map_func_type
        assert mapper.nullvalue == self.nullval


    def test_map(self):

        a, b, c = self.coordnames
        self.mapper.map(self.vtu_out, a, b, c)
        self.mapper.write_vtp(self.vtp_out)

        # check variable names
        assert isinstance(self.mapper.nc["coord_names"], dict)
        assert self.mapper.nc["coord_names"]["lat_name"] == self.coordnames[0]
        assert self.mapper.nc["coord_names"]["lon_name"] == self.coordnames[1]
        assert self.mapper.nc["coord_names"]["time_name"] == self.coordnames[2]

        assert self.mapper.nc["data_names"] == self.nc_data_var

        # check coordinate data
        self.mapper.nc["coords"]["lat"] = self.lat_dat
        self.mapper.nc["coords"]["lon"] = self.lon_dat
        self.mapper.nc["coords"]["time"] = self.time_dat

        # check extracted data
        assert_array_almost_equal(
            self.mapper.nc["data"][0][1],
            self.nc_data[0][1]
        )

        assert_array_almost_equal(
            self.mapper.nc["data"][1][1],
            self.nc_data[1][1]
        )

        # check created vtp and vtu
        isinstance(self.mapper.vtp, vtkPolyData)
        isinstance(self.mapper.out_vtu, vtkUnstructuredGrid)

        # check written files
        isfile(self.vtu_out)
        isfile(self.vtp_out)
