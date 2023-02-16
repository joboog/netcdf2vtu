"""A Python package to interpolate data from a netcdf file onto
a VTU file.

**netcdf2vtu** provides a the *Mapper* class to read in a netcdf input
file, a destination VTU file and to set some specifics for data
interpolation.
"""

from netCDF4 import Dataset
import numpy as np
from pyproj import Transformer
import vtk
from vtk.util import numpy_support


class Mapper(object):
    """
    A class for reading and interpolating data from a netcdf file onto
    a VTU file.

    Parameters
    ----------
    nc_path : str
        The path to the netcdf input file.
    vtu_path : str
        The path to the destination VTU file.
    data_var_names : list of str
        A list of variable names to read from the netcdf file.
    map_func_type : int
        The type of mapping function to use for interpolation.
    nc_crs : str
        The coordinate reference system of the netcdf file.
    vtu_crs : str
        The coordinate reference system of the VTU file.
    nullvalue : float
        The value to use for missing data points during interpolation.
    kwargs_nc : dict
        Additional keyword arguments to pass to the netcdf4 library when
         reading the netcdf file.

    """

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
        """
        Initializes the Mapper object with the provided inputs.
        """

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
        """
        Returns a numpy.ndarray variable data of the netcdf dataset.
        """

        self.nc["dataset"].variables


    def set_nc_coord_names(self,
                            lat_name,
                            lon_name,
                            time_name = None):
        """
        Sets the names of the latitude, longitude, and time variables
        to be read from the dataset in `self.nc["dataset"]`.
        """

        self.nc.update({"coord_names" : {"lat_name" : lat_name,
                                        "lon_name" : lon_name,
                                        "time_name" : time_name}})


    def set_nc_data_names(self, var_names):
        """
        Sets the names of the variables to be read from the netcdf
        dataset in `self.nc["dataset"]`.

        Parameter
        ---------
        var_names : list of str
            A list of variable names to read from the netcdf file.
        """

        self.nc["data_names"] = var_names


    def read_nc_coords(self):
        """
        Reads the latitude, longitude, and time coordinates from the
        netcdf dataset in `self.nc["dataset"]`.
        """

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
        """
        Reads the data variables specified in `self.nc["data_names"]`
        from the netcdf dataset.
        """

        nc_data = get_nc_data(self.nc["dataset"],
                                  self.nc["data_names"])
        self.nc.update({"data" : nc_data})


    def interpolate(self):
        """
        Interpolates the data from the netcdf file onto the VTU file.
        It creates a vtkPolyData instance from `self.nc["dataset"]` and
        interpolates this on a copy of an instance of
        vtkUnstructuredGrid created from the input VTU file.
        """

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
        """
        Writes the updated instance of vtkUnstructuredGrid to the
        specified path as VTU file.

        Parameter
        ---------
        path : str
            The path of the VTU file to be written.
        """

        write_vtu(self.out_vtu, path)
        print("New VTU mesh written to disk.")


    def write_vtp(self, path):
        """
        Writes the vtkPolyData instance created to the
        specified path as VTP file.

        Parameter
        ---------
        path : str
            The path of the VTP file to be written.
        """

        write_vtp(self.vtp, path)


    def map(self,
            out_path,
            lat_name,
            lon_name,
            time_name = None):
        """
        A wrapper for simple interpolation.
        """

        self.set_nc_coord_names(lat_name,
                                lon_name,
                                time_name)
        self.read_nc_coords()
        self.read_nc_data()
        self.interpolate()
        self.write_out_vtu(out_path)



def get_nc_data(nc, data_var_names):
    """
    Extract variable data out of a netcdf4 dataset.

    Parameters
    ----------
    nc : netCDF4.Dataset
        The netcdf dataset to extract data from.
    data_var_names : list of str
        A list of variable names to extract from nc.

    Returns
    -------
    list
        A list of tuples, where each tuple contains a variable name and the corresponding data as a numpy array.
    """

    vars_data = []
    for i, var_name in enumerate(data_var_names):
        t = (var_name, nc.variables[var_name][:].filled())
        vars_data.append(t)

    return(vars_data)


def create_vtp(lon_dat, lat_dat, nc_crs, vtu_crs):
    """
    Initializes a vtkPolyData object and sets its points and cells based
     on input latitude and longitude coordinates.
    The point coordinates are transformed to the specified destination
    coordinate reference system.

    initialize src_poly as vtkPolyData
    set points and cells for src_poly (vtkPolyData object where netcdf
    data goes into)
    cells are vertice cells based on points
    point coordinates are transformed

    Parameters
    ----------
    lon_dat : numpy.ndarray
        A numpy array of longitude coordinates.
    lat_dat : numpy.ndarray
        A numpy array of latitude coordinates.
    nc_crs : str
        The coordinate reference system of the input latitude and
        longitude coordinates.
    vtu_crs : str
        The destination coordinate reference system to transform the
        point coordinates to.

    Returns
    -------
    vtk.vtkPolyData
        A vtkPolyData object with points and cells set based on input
        coordinates.
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
    Adds extracted variable data from a netcdf file to a vtkPolyData
    object.

    Parameters
    ----------
    vtp : vtk.vtkPolyData
        The vtkPolyData object to add the netcdf data to.
    nc_data : list
        A list of tuples, where each tuple contains a variable name and its corresponding data as a numpy array.
    time : numpy.ndarray, optional
        An array of time values for the netcdf data. If not provided,
        the time value will be set to "FALSE".
    """

    if time is None:
        time = "F"

    for j,val in enumerate(nc_data):
        for i in range(len(time)):

            if type(time) is str:
                arr_name = val[0]
            else:
                arr_name = val[0] + "_%s" % str(int(time[i]))

            new_point_arr_vtk = numpy_support.numpy_to_vtk(val[1][i].ravel())
            new_point_arr_vtk.SetName(arr_name)
            vtp.GetPointData().AddArray(new_point_arr_vtk)

    vtp.Modified()



def write_vtp(vtp, outputfile_name):
    """
    Writes a vtkPolyData object to a file.

    Parameters
    ----------
    vtp : vtk.vtkPolyData
        The vtkPolyData object to write to a file.
    outputfile_name : str
        The path of the file to write the vtkPolyData object to.
    """

    write_ouput = vtk.vtkXMLPolyDataWriter()
    write_ouput.SetInputData(vtp)
    write_ouput.SetFileName(outputfile_name)
    write_ouput.Write()


def read_vtu(in_filepath):
    """
    Reads an OGS mesh file and returns it as a vtkUnstructuredGrid object.

    Parameters
    ----------
    in_filepath : str
        The path of the VTU file to read.

    Returns
    -------
    vtk.vtkUnstructuredGrid
        The VTU file as a vtkUnstructuredGrid instance.
    """

    dst = vtk.vtkXMLUnstructuredGridReader()
    dst.SetFileName(in_filepath)
    dst.Update()
    return(dst.GetOutput())


def interpolate_vtp_data_on_vtu(vtp, vtu, map_func_type, nullvalue):
    """
    Maps the data of a vtkPolyData object onto a vtkUnstructuredGrid
    object using a specified interpolation algorithm.

    Parameters
    ----------
    vtp : vtk.vtkPolyData
        The object that contains the data to be interpolated.
    vtu : vtk.vtkUnstructuredGrid
        The object that the data will be interpolated onto.
    map_func_type : int
        The type of mapping function to use for interpolation.
    nullvalue : float
        The value to use for missing data points during interpolation.

    Returns
    -------
    vtk.vtkUnstructuredGrid
        An object with the interpolated data.
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
    Writes a vtkUnstructuredGrid  object to a file.

    Parameters
    ----------
    vtp : vtk.vtkUnstructuredGrid
        The vtkUnstructuredGrid object to write to a file.
    outputfile_name : str
        The path of the VTU file to be written.
    """
    write_output = vtk.vtkXMLUnstructuredGridWriter()
    write_output.SetFileName(path)
    write_output.SetInputData(vtu)
    write_output.Write()


def set_map_function(map_func_type):
    """
    Returns the specified interpolation/mapping algorithm.

    Parameters
    ----------
    map_func_type : int
        The type of mapping function to use for interpolation.
        1: voronoi_kernel
        2: gaussian_kernel
        3: shepard_kernel

    Returns
    -------
    vtk.vtkInterpolationKernel
        An object of the specified interpolation/mapping algorithm.
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
    Gaussian filter of point cloud defined by
    radius around a point.

    Returns
    -------
    vtk.vtkGaussianKernel
        The Gaussian filter interpolation kernel.
    """
    int_kernel = vtk.vtkGaussianKernel()
    int_kernel.SetSharpness(2)
    int_kernel.SetRadius(4000)
    #TODO: add kwargs to set the sharpness and radius
    return(int_kernel)


def voronoi_kernel():
    """
    Preferred for categorical data.
    Takes value of the closest point.

    Returns
    -------
    vtk.vtkVoronoiKernel
        The voronoi interpolation kernel.
    """
    int_kernel = vtk.vtkVoronoiKernel()
    return(int_kernel)


def shepard_kernel():
    """
    Interpolation of a point cloud object defined by
    power function of the radius around a point.

    Returns
    -------
    vtk.vtkShepardKernel
        The Shepard interpolation kernel.
    """
    int_kernel = vtk.vtkShepardKernel()
    int_kernel.SetPowerParameter(2)
    int_kernel.SetRadius(4000)
    return(int_kernel)