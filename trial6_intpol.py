#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 12:46:21 2019

@author: boog
"""

import vtk
#import numpy as np
from vtk.util import numpy_support
from gstools import SRF, Gaussian  

# =============================================================================
# import src data
# =============================================================================
src = vtk.vtkXMLStructuredGridReader()
src.SetFileName("/home/boog/ufz/08_hiwi_envinf/07_netcdf2vtk" \
                   "/01_paraview_trial/01_interpol_accuracy" \
                   "/01_meshes/src_4000x4000.vts")
src.Update()
src = src.GetOutput()


# =============================================================================
# generate dummy src data
# =============================================================================

# structured field with a size of 100x100 and a grid-size of 1x1
x = y = range(100)
#x = y = np.arange(0,400000,4000)
model = Gaussian(dim=2, var=1, len_scale=10, nugget=0)
srf = SRF(model)
srf.structured([x, y])

# tranfrom to vtk format and add as array data to output
# use .flatten() to convert 2d array to 1d
new_point_arr_vtk = numpy_support.numpy_to_vtk(srf.field.flatten())
new_point_arr_vtk.SetName("rndgaussian")
src.GetPointData().AddArray(new_point_arr_vtk)


# =============================================================================
# import dst mesh
# =============================================================================
dst = vtk.vtkXMLUnstructuredGridReader()
dst.SetFileName("/home/boog/ufz/08_hiwi_envinf/07_netcdf2vtk" \
                   "/01_paraview_trial/01_interpol_accuracy" \
                   "/01_meshes/dst_40000x40000.vtu")
dst.Update()
dst = dst.GetOutput()


# =============================================================================
# map data
# =============================================================================

# def interpolation kernel
gaussian_kernel = vtk.vtkGaussianKernel()
gaussian_kernel.SetSharpness(2)
gaussian_kernel.SetRadius(4000)

# def interpolion function 
interpolator = vtk.vtkPointInterpolator()
interpolator.SetInputData(dst)
interpolator.SetSourceData(src)
interpolator.SetKernel(gaussian_kernel)
interpolator.SetNullValue(-9999)
interpolator.Update()

int_out = interpolator.GetOutput()        
        

# =============================================================================
# write output mesh
# =============================================================================
write_ouput = vtk.vtkXMLUnstructuredGridWriter()
write_ouput.SetFileName('trial6.vtu')
write_ouput.SetInputData(int_out)
write_ouput.Write()