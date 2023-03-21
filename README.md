# NetCDF2VTU

[![Pipeline Status](https://gitlab.com/joboog/netcdf2vtu/badges/master/pipeline.svg)](https://gitlab.com/joboog/netcdf2vtu/-/commits/master)
[![Test Coverage](https://gitlab.com/joboog/netcdf2vtu/badges/master/coverage.svg)](https://gitlab.com/joboog/netcdf2vtu/-/jobs)
[![Documentation Status](https://readthedocs.org/projects/netcdf2vtu/badge/?version=latest)](https://netcdf2vtu.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/netcdf2vtu.svg)](https://badge.fury.io/py/netcdf2vtu)

Netcdf2vtu is a Python package including a command line interface (CLI)
 to interpolate data from netCDF on VTU files.

## Setup

Ideally, you have created a virtual python environment. If not look
[here](https://packaging.python.org/en/latest/tutorials/installing-packages/#creating-virtual-environments).

A stable release is available on https://pypi.org so you can install it via
```
$ python3 -m pip install netcdf2vtu
```

To install the development version, you have to actually install it from a clone of this repo.
```
$ git clone https://gitlab.com/joboog/netcdf2vtu.git
```

Go into the repo and install the dependencies:
```
$ cd netcdf2vtu
$ pip install -r requirements.txt
```

Then install the package:
```
$ pip install .
```

You can run the tests with:
```
$ pytest .
```

## Usage

Netcdf2vtu can be imported as a Python module or used from the command
 line.


### As imported module

A brief example of how to map data from the NETCDF4 file `data/ex1_3.nc`
on the destination VTU `data/ogs.vtu` using the `Mapper` class is
 presented below.
Just define the required input data to the `Mapper` class.

```
from netcdf2vtu.netcdf2vtu import Mapper

# setup ---------------------------------------------------------
nc_path = "data/ex1_3.nc" # path of input netcdf file
vtu_path = "data/ogs.vtu" # path of vtu file
vtu_new_path = "ex3_new.vtu" # path of updated vtu file
data_var_names = ["SWC_L01", "SWC_L02"] # names of the netcdf data \
                                        # variables
map_func_type = 1 # def mapping func type 1: Voronoi, 2:Gaussian,
                  # 3:Shepard
nc_crs = "EPSG:4326"   # coordinate system in netcdf file
vtu_crs = "EPSG:5684"   # coordinate systetm in vtu file
nullvalue = -9999
```

Then create a mapper object.

```
# for convenience, all the above can be within two statements
mapper = Mapper(nc_path,
                vtu_path,
                data_var_names,
                map_func_type,
                nc_crs,
                vtu_crs,
                nullvalue)
```

And start the interpolation:
```
mapper.map(out_path = vtu_new_path,
           lat_name = "lat",
           lon_name = "lon",
           time_name = "time")
```

The outputted file `ex3_new.vtu` is a copy of the VTU file
`data/ogs.vtu` including the interpolated data.


### Command line interface

From the command line you can do the same as above with:

```
$ netcdf2vtu -o ex3_new.vtu --time time data/ex1_3.nc data/ogs.vtu EPSG:4326 EPSG:5684 SWC_L01 SWC_L02
```
Note that the arguments to define the name of the output VTU file, the `map_func_type` and names for coordinates have default values and may not need to be explicitly set in the call as shown above.
Use the help function to see the specifics of the `netcdf2vtu` command:

```
$ netcdf2vtu -h

usage: netcdf2vtu [-h] [-o OUT_VTU] [-m MAP_TYPE] [-n NULLVALUE] [--lat LAT] [--lon LON] [--time TIME] in_nc in_vtu nc_crs vtu_crs var_names [var_names ...]

Interpolates data from netCDF4 files to VTU files.

positional arguments:
  in_nc                 Path to the input netCDF file.
  in_vtu                Path to the input VTU file.
  nc_crs                CRS of input netcdf file
  vtu_crs               CRS of input vtu file.
  var_names             Names of the data variables in the input netcdf file.

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_VTU, --out_vtu OUT_VTU
                        Path to the output VTU file; default is "out.vtu".
  -m MAP_TYPE, --map_type MAP_TYPE
                        Type of interpolation function; default is 1.
  -n NULLVALUE, --nullvalue NULLVALUE
                        Nullvalue; default is -9999.
  --lat LAT             Name of the latitude variable; default is "lat".
  --lon LON             Name of the longitude variable; default is "lon".
  --time TIME           Name of the time variable if present.
```

## Contribution

If there are any problems please raise an issue [here](https://gitlab.com/joboog/netcdf2vtu/-/issues).

Please look at [CONTRIBUTING.md](https://gitlab.com/joboog/netcdf2vtu/-/blob/master/CONTRIBUTING.md).

