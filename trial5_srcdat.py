#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 12:46:21 2019

@author: boog
"""

import vtk
import numpy as np
from vtk.util import numpy_support
from gstools import SRF, Gaussian  

# import vts
reader = vtk.vtkXMLStructuredGridReader()
reader.SetFileName("/home/boog/ufz/08_hiwi_envinf/07_netcdf2vtk" \
                   "/01_paraview_trial/01_interpol_accuracy" \
                   "/01_meshes/src_4000x4000.vts")
reader.Update()
output = reader.GetOutput()


# add point data array
# structured field with a size of 100x100 and a grid-size of 1x1
x = y = range(100)
#x = y = np.arange(0,400000,4000)
model = Gaussian(dim=2, var=1, len_scale=10, nugget=0)
srf = SRF(model)
srf.structured([x, y])

# tranfrom to vtk format and add as array data
# use .flatten() to convert 2d array to 1d
new_point_arr_vtk = numpy_support.numpy_to_vtk(srf.field.flatten())
new_point_arr_vtk.SetName("rndgaussian")
output.GetPointData().AddArray(new_point_arr_vtk)

# write as new vts
write_ouput = vtk.vtkXMLStructuredGridWriter()
write_ouput.SetFileName('trial5_io.vts')
write_ouput.SetInputData(output)
write_ouput.Write()