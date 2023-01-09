#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 12:46:21 2019

@author: boog
"""

import vtk
import numpy as np
from vtk.util import numpy_support

# import vts
reader = vtk.vtkXMLStructuredGridReader()
reader.SetFileName("/home/boog/ufz/08_hiwi_envinf/07_netcdf2vtk" \
                   "/01_paraview_trial/01_interpol_accuracy" \
                   "/01_meshes/src_4000x4000.vts")
reader.Update()

output = reader.GetOutput()

# add point data array
new_point_arr0 = np.empty([100,100])
x = np.arange(0,400000,4000, dtype='Float64')

for i in xrange(0,len(new_point_arr0)):
    new_point_arr0[i] = x
        

# tranfrom to vtk format and add as array data
# use .flatten() to convert 2d array to 1d
new_point_arr_vtk = numpy_support.numpy_to_vtk(new_point_arr0.flatten())
output.GetPointData().AddArray(new_point_arr_vtk)

# write as new vts
write_ouput = vtk.vtkXMLStructuredGridWriter()
write_ouput.SetFileName('trial3_io.vts')
write_ouput.SetInputData(output)
write_ouput.Write()