# NetCDF2VTU

Netcdf2vtu is a Python package to interpolate data from netCDF on VTU
 files.

## Setup

The package is still under development and not yet available on as
package on https://pypi.org .
You have to actually install from a clone of the repo.

Ideally, you have created a virtual python environment. If not look
[here](https://packaging.python.org/en/latest/tutorials/installing-packages/#creating-virtual-environments).

Then clone the repo:
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
$ pytest
```

## Usage

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
mapper.map(out_path = "ex3_mapper2_new_vtu.vtu",
           lat_name = "lat",
           lon_name = "lon",
           time_name = "time")
```

The outputted file `ex3_mapper2_new_vtu.vtu` is a copy of the VTU file
including the interpolated data.

## Contribution

If there are any problems please raise an issue [here](https://gitlab.com/joboog/netcdf2vtu/-/issues).)

Please look at [CONTRIBUTING.md](./CONTRIBUTING.md).

