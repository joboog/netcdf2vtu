"""Command line interface for netcdf2vtu."""

import argparse

from .netcdf2vtu import Mapper

def main():

    # Define parser
    parser = argparse.ArgumentParser(
        prog="netcdf2vtu",
        description='Interpolates data from  netCDF4 files to VTU files.')

    parser.add_argument('-o', '--out_vtu', type=str, default="out.vtu",
                        help='Path to the output VTU file; default is '\
                            '"out.vtu".')
    parser.add_argument('-m', '--map_type', type=int, default=1,
                        help='Type of interpolation function; default is'\
                            ' 1.')
    parser.add_argument('-n', '--nullvalue', type=int, default=-9999,
                        help='Nullvalue; default is -9999.')
    parser.add_argument('--lat', type=str, default="lat",
                        help='Name of the latitude variable; default' \
                            ' is "lat".')
    parser.add_argument('--lon', type=str, default="lon",
                        help='Name of the longitude variable; default' \
                            ' is "lon".')
    parser.add_argument('--time', type=str,
                        help='Name of the time variable if present.')
    parser.add_argument('in_nc', type=str,
                        help='Path to the input netCDF file.')
    parser.add_argument('in_vtu', type=str,
                        help='Path to the input VTU file.')
    parser.add_argument('nc_crs', type=str,
                        help='CRS of input netcdf file')
    parser.add_argument('vtu_crs', type=str,
                        help='CRS of input vtu file.')
    parser.add_argument('var_names', nargs="+", type= str,
                        help='Names of the data variables in the input netcdf file.')

    args = parser.parse_args()

    # Assign args
    in_nc = args.in_nc
    in_vtu = args.in_vtu
    out_vtu = args.out_vtu
    var_names = args.var_names
    map_type = args.map_type
    nc_crs = args.nc_crs
    vtu_crs = args.vtu_crs
    nullval = args.nullvalue
    lat_name = args.lat
    lon_name = args.lon
    time_name = args.time

    # Do the interpolation
    mapper = Mapper(in_nc,
                    in_vtu,
                    var_names,
                    map_type,
                    nc_crs,
                    vtu_crs,
                    nullval)

    mapper.map(out_path = out_vtu,
                lat_name = lat_name,
                lon_name = lon_name,
                time_name = time_name)

if __name__ == "__main__":
    main()