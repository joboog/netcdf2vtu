#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: joboog
"""


import numpy as np
import vtk
from vtk.util import numpy_support
from pyproj import Transformer

from netcdf2vtk.mapping import set_map_function




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