#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: joboog
"""

from netCDF4 import Dataset
import numpy as np
from pyproj import Transformer
import vtk
from vtk.util import numpy_support


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
        dst_vtu = read_vtu(self.in_vtu["path"])
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

        nc_data = get_nc_data(self.nc["dataset"],
                                  self.nc["data_names"])
        self.nc.update({"data" : nc_data})


    def interpolate(self):

        self.vtp = create_vtp(
                    self.nc["coords"]["lon"],
                    self.nc["coords"]["lat"],
                    self.nc["crs"],
                    self.in_vtu["crs"])

        nc_data_to_vtp(self.vtp,
                           self.nc["data"],
                           self.nc["coords"]["time"])

        print("NetCDF converted to VTP.")


        self.out_vtu = interpolate_vtp_data_on_vtu(
                        self.vtp,
                        self.in_vtu["vtu"],
                        self.map_func_type,
                        self.nullvalue)

        print("Data from VTP interpolated to VTU.")


    def write_out_vtu(self, path):
        # output updated ogs-mesh
        write_vtu(self.out_vtu, path)
        print("New VTU mesh written to disk.")


    def write_vtp(self, path):

        write_vtp(self.vtp, path)


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



def get_nc_data(nc, data_var_names):
    """
    extract variable data out of imported netcdf (src_nc)
    """
    var_data = [None] * len(data_var_names)
    vars = np.array([data_var_names, var_data]).T
    for i in range(len(vars)):
        vars[i][1] = nc.variables[(vars[i][0])][:].filled()
    return(vars)


def create_vtp(lon_dat, lat_dat, nc_crs, vtu_crs):
    """
    initialize src_poly as vtkPolyData
    set points and cells for src_poly (vtkPolyData object where netcdf
    data goes into)
    cells are vertice cells based on points
    point coordinates are transformed

    returnes src_poly
    """
    # check dims and def point array
    if len(lon_dat.shape) == 1 and len(lat_dat.shape) == 1:

        xv, yv, zv = np.meshgrid(lon_dat, lat_dat, np.array([0]))
        xv = xv.ravel()
        yv = yv.ravel()
        zv = zv.ravel()

    elif len(lon_dat.shape) == 2 and len(lat_dat.shape) == 2:

        xv = lon_dat.ravel()
        yv = lat_dat.ravel()
        zv = np.repeat(0, len(xv))

    else:
        raise ValueError("Shape of coordinate arrays not as expected!")

    points = np.array([xv, yv, zv]).transpose()

    # def vtkPoints/cells
    nc_points = vtk.vtkPoints()
    nc_cells = vtk.vtkCellArray()

    # transform points coordiantes
    transf = Transformer.from_crs(
                            nc_crs,
                            vtu_crs,
                            always_xy=True)

    for i in range(len(points)):
        p = transf.transform(points[i][0], points[i][1])
        ind = nc_points.InsertNextPoint(p[0], p[1], 0)
        nc_cells.InsertNextCell(1)
        nc_cells.InsertCellPoint(ind)

    # define vtkPolydata obj
    vtp = vtk.vtkPolyData()
    vtp.SetPoints(nc_points)
    vtp.SetVerts(nc_cells)
    return(vtp)


def nc_data_to_vtp(vtp, nc_data, time = None):
    """
    adds extracted data arrays from imported netcdf (src_nc) to
    vtkPolyData obj (src_poly)
    """
    if time is None:
        time = "F"

    for j in range(len(nc_data)):
        for i in range(len(time)):

            if type(time) is str:
                arr_name = nc_data[j][0]
            else:
                arr_name = nc_data[j][0] + "_%s" % str(int(time[i]))

            new_point_arr_vtk = numpy_support.numpy_to_vtk(nc_data[j][1][i].ravel())
            new_point_arr_vtk.SetName(arr_name)
            vtp.GetPointData().AddArray(new_point_arr_vtk)

    vtp.Modified()



def write_vtp(vtp, outputfile_name):
    """
    writes imported netcdf as vtkPolyData obj
    """
    write_ouput = vtk.vtkXMLPolyDataWriter()
    write_ouput.SetInputData(vtp)
    write_ouput.SetFileName(outputfile_name)
    write_ouput.Write()


def read_vtu(in_filepath):
    """
    reads ogs mesh file where the netcdf data will be mapped on
    """
    dst = vtk.vtkXMLUnstructuredGridReader()
    dst.SetFileName(in_filepath)
    dst.Update()
    return(dst.GetOutput())


def interpolate_vtp_data_on_vtu(vtp, vtu, map_func_type, nullvalue):
    """
    maps the data of src_poly on imported ogs-mesh (ogs_vtu) using the
    interpolation algorithm defined in set_int_kernel()
    """
    interpolator = vtk.vtkPointInterpolator()
    interpolator.SetInputData(vtu)
    interpolator.SetSourceData(vtp)
    # def interpolation algorithm
    interpolator.SetKernel(set_map_function(map_func_type))
    # def value if interpolation does not work
    interpolator.SetNullValue(nullvalue)
    interpolator.Update()
    return(interpolator.GetOutput())



def write_vtu(vtu, path):
    """
    write the ogs-mesh including the newly mapped data
    """
    write_output = vtk.vtkXMLUnstructuredGridWriter()
    write_output.SetFileName(path)
    write_output.SetInputData(vtu)
    write_output.Write()


def set_map_function(map_func_type):
    """
    defines the mapping/interpolation algorithm
    """
    switcher = {
        1: voronoi_kernel,
        2: gaussian_kernel,
        3: shepard_kernel
    }
    # Get the function from switcher dictionary
    func = switcher.get(map_func_type, lambda: "Invalid func_type")
    return(func())


def gaussian_kernel():
    """
    gaussian filter of point cloud in src_poly defined by radius
    around point in ogs_vtu
    """
    int_kernel = vtk.vtkGaussianKernel()
    int_kernel.SetSharpness(2)
    int_kernel.SetRadius(4000)
    return(int_kernel)


def voronoi_kernel():
    """
    preferred for categorial data
    takes value of closest points in src_poly
    """
    int_kernel = vtk.vtkVoronoiKernel()
    return(int_kernel)


def shepard_kernel():
    """
    interpolation of point cloud in src_poly defined by power of radius
    around point in ogs_vtu
    """
    int_kernel = vtk.vtkShepardKernel()
    int_kernel.SetPowerParameter(2)
    int_kernel.SetRadius(4000)
    return(int_kernel)