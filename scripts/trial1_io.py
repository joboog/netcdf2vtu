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

# add point data array
output = reader.GetOutput()
points = output.GetPointData() 

new_point_arr = np.arange(0,400000,400, dtype='Float64')
new_point_arr_vtk = numpy_support.numpy_to_vtk(new_point_arr)

output.GetPointData().AddArray(new_point_arr_vtk)

# write as new vts
write_ouput = vtk.vtkXMLStructuredGridWriter()
write_ouput.SetFileName('trial_output.vts')
write_ouput.SetInputData(output)
write_ouput.Write()