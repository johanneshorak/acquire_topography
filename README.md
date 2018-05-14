# acquire_topography
 script name: act.py
 author	    : Johannes Horak
 
A script that generates an upscaled topography of a region from a high resolution DEM

# description:
act.py is a script that may be used to generate a digital elevation model (DEM) in lambert conformal projection for downscaling purposes for a specified region of the world with a uniform horizontal grid spacing. If so desired topographic smoothing may be applied by using the additional option --smooth. See filter_topo for more details. Currently only DEMs with a minimum grid spacing of 1' may be	generated.

All parameters required are entered via the command line.

# some examples:
    python act.py --lat1=45 --lat2=48 --lat0=48 --lonc=10 --nx=205 --ny=150 --dx 4000 --dy 4000 --name europe
generates a DEM for the European Alps.

    python act.py --lat1=-42 --lat2=-46 --lat0=-43.6 --lonc=170 --nx=205 --ny=225 --dx 4000 --dy 4000 --name southern_nz --smooth
generates a smoothed DEM of the South Island of New Zealand.

# requirements:
The ETOPO1 DEM is required to be placed in the [script-dir]/data directory in the NETCDF format. It may be downloaded from https://www.ngdc.noaa.gov/mgg/global/global.html. Adjust the filename of the ETOPO1-Ice DEM accordingly by changing	the value of the variable ETOPO1I_FN.

# packages required:
xarray, numpy, scipy, matplotlib (for debug or preview plots), salem, bunch, getopt
