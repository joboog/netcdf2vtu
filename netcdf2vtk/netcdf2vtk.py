#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: joboog
"""

from netCDF4 import Dataset
import numpy as np
import vtk
from vtk.util import numpy_support
from pyproj import Transformer

from netcdf2vtk.mapping import set_map_function


def get_src_nc_data(src_nc, data_var_names):
    """
    extract variable data out of imported netcdf (src_nc)
    """
    var_data = [None] * len(data_var_names)
    src_vars = np.array([data_var_names, var_data]).T
    for i in range(len(src_vars)):
        src_vars[i][1] = src_nc.variables[(src_vars[i][0])][:].filled()
    return(src_vars)


def init_src_poly(lon_dat, lat_dat, src_crs, dst_crs):
    """
    initialize src_poly as vtkPolyData
    set points and cells for src_poly (vtkPolyData object where netcdf
    data goes into)
    cells are vertice cells based on points
    point coordinates are transformed

    returnes src_poly
    """
    # def point array
    points = np.array([lon_dat.flatten(),
                       lat_dat.flatten(),
                       np.repeat(0, len(lat_dat.flatten()))
                    ]).transpose()

    # def vtkPoints/cells
    src_points = vtk.vtkPoints()
    src_cells = vtk.vtkCellArray()

    # transform points coordiantes
    transf = Transformer.from_crs(
                            src_crs,
                            dst_crs,
                            always_xy=True)

    for i in range(len(points)):
        p = transf.transform(points[i][0], points[i][1])
        ind = src_points.InsertNextPoint(p[0], p[1], 0)
        src_cells.InsertNextCell(1)
        src_cells.InsertCellPoint(ind)

    # define vtkPolydata obj
    src_poly = vtk.vtkPolyData()
    src_poly.SetPoints(src_points)
    src_poly.SetVerts(src_cells)
    return(src_poly)


def add_nc_data_to_src_poly(src_poly, time, src_vars):
    """
    adds extracted data arrays from imported netcdf (src_nc) to
    vtkPolyData obj (src_poly)
    """
    for j in range(len(src_vars)):
        for i in range(len(time)):
            arr_name = src_vars[j][0] + "_%s" % str(int(time[i]))
            new_point_arr_vtk = numpy_support.numpy_to_vtk(src_vars[j][1][i].flatten())
            new_point_arr_vtk.SetName(arr_name)
            src_poly.GetPointData().AddArray(new_point_arr_vtk)

    src_poly.Modified()



def write_src_poly(src_obj, outputfile_name):
    """
    writes imported netcdf as vtkPolyData obj
    """
    write_ouput = vtk.vtkXMLPolyDataWriter()
    write_ouput.SetInputData(src_obj)
    write_ouput.SetFileName(outputfile_name)
    write_ouput.Write()


def read_ogs_vtu(in_filepath):
    """
    reads ogs mesh file where the netcdf data will be mapped on
    """
    dst = vtk.vtkXMLUnstructuredGridReader()
    dst.SetFileName(in_filepath)
    dst.Update()
    return(dst.GetOutput())


def map_data_on_ogs_vtu(src_poly, ogs_vtu, map_func_type, nullvalue):
    """
    maps the data of src_poly on imported ogs-mesh (ogs_vtu) using the
    interpolation algorithm defined in set_int_kernel()
    """
    interpolator = vtk.vtkPointInterpolator()
    interpolator.SetInputData(ogs_vtu)
    interpolator.SetSourceData(src_poly)
    # def interpolation algorithm
    interpolator.SetKernel(set_map_function(map_func_type))
    # def value if interpolation does not work
    interpolator.SetNullValue(nullvalue)
    interpolator.Update()
    return(interpolator.GetOutput())



def write_mapped_ogs_vtu(out_vtu, out_filename):
    """
    write the ogs-mesh including the newly mapped data
    """
    write_output = vtk.vtkXMLUnstructuredGridWriter()
    write_output.SetFileName(out_filename)
    write_output.SetInputData(out_vtu)
    write_output.Write()