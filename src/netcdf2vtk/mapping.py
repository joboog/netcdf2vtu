#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: joboog
"""

import vtk

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