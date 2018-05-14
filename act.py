# ----------------------------------------------------------------------
# script name: act.py
# author	 : Johannes Horak
#
# description:
#	act.py is a script that may be used to generate a digital elevation
#	model (DEM) in lambert conformal projection for downscaling purposes
#	for a specified region of the world with a uniform horizontal
#	grid spacing.
#
#	If so desired topographic smoothing may be applied by using the
#	additional option --smooth. See filter_topo for more details.
#
#	Currently only DEMs with a minimal grid spacing of 1' may be
#	generated.
#	
#	All parameters required are entered via the command line.
#
#   some examples:
#		python act.py --lat1=45 --lat2=48 --lat0=48 --lonc=10 --nx=205 --ny=150 --dx 4000 --dy 4000 --name europe
#	generates a DEM of the European Alps and some of the surrounding
#	regions.
#
#		python act.py --lat1=-42 --lat2=-46 --lat0=-43.6 --lonc=170 --nx=205 --ny=225 --dx 4000 --dy 4000 --name southern_nz --smooth
#	generates a smoothed DEM of the South Island of New Zealand
#	including the surrounding ocean.
#
# requirements:
#	The ETOPO1 DEM is required to be placed in the [script-dir]/data
#	directory in the NETCDF format. It may be downloaded from
#	https://www.ngdc.noaa.gov/mgg/global/global.html
#	Adjust the filename of the ETOPO1-Ice DEM accordingly by chaning
#	the value of the variable ETOPO1I_FN.
#
# ----------------------------------------------------------------------
import math   	  as math
import xarray     as xa
import numpy  	  as np
import pandas  	  as pa
import scipy   	  as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import salem
import scipy.interpolate as interpolate
from bunch     import Bunch
from math      import radians, cos, sin, asin, sqrt
import warnings
import getopt
import glob, os, sys

warnings.simplefilter(action='ignore', category=FutureWarning)

R 			= 6371 						# Radius of earth in kilometers.
ETOPO1I_FN	= 'ETOPO1_Ice_g_gmt4.grd'	# change this if your ETOPO1-Ice DEM file is named differently

# ----------------------------------------------------------------------
# d2r(angle)
#
# short: converts an angle from degrees to radians
#
# arguments:
#	angle	... angle in degree
#
# returns:
#	angle converted to radians
#
# ----------------------------------------------------------------------
def d2r(angle):
	return angle*np.pi/180.0

# ----------------------------------------------------------------------
# haversine(lon1, lat1, lon2, lat2)
#
# short: calculates the distance in km bzw. two points on the surface
# 		 of the earth.
#
# arguments:
#	lon1	... longitude of point 1 on earth in degrees
#	lat1	... latitude of point 1 on earth in degrees
#	lon2	... longitude of point 2 on earth in degrees
#	lat2	... latitude of point 2 on earth in degrees
#
# returns:
# 	the distance between point 1 and 2 in kilometers
# ----------------------------------------------------------------------
def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    return c * R
    
# ----------------------------------------------------------------------
# filter_topo(ds)
#
# short: applies a 3x3 moving windows smoothing to a DEM
#
# arguments:		
#	ds		... a dataset containing  the high resolution DEM
#
# returns:
#	ds_new	... the filtered version of the input topography.
#
# description:
#	the function applies a 3x3 moving window smoothing algorithm
#	that essentially removes peaks that are represented only by
#	a single grid cell.
#
#	reference:  
#	   GUO, Yong-run and SUE CHEN, Terrain and land use for the 
#	   fifth-generation Penn State/NCAR Mesoscale Modeling System (MM5):
#	   Program TERRAIN. 1994. 
#
# ----------------------------------------------------------------------
def filter_topo(ds):
	smcf=0.5
	
	nlon_lr = len(ds.west_east)
	nlat_lr = len(ds.south_north)

	htmp = np.empty(nlon_lr*nlat_lr)
	htmp = htmp.reshape(nlat_lr,nlon_lr)
	
	hsmth = np.empty(nlon_lr*nlat_lr)
	hsmth = hsmth.reshape(nlat_lr,nlon_lr)
	
	for i in range(0,nlat_lr):
		for j in range(0,nlon_lr):
			if (i > 0 and j> 0 and i<(nlat_lr-1) and j<(nlon_lr-1)):
				#print str(i).zfill(4),",",str(j).zfill(4)
				htmp[i,j]=ds.HGT_M[i,j]+smcf*(0.5*(ds.HGT_M[i,j+1]+ds.HGT_M[i,j-1])-ds.HGT_M[i,j])
			else:
				htmp[i,j]=ds.HGT_M[i,j]
		
	for i in range(0,nlat_lr):
		for j in range(0,nlon_lr):
			if (i > 0 and j> 0 and i<(nlat_lr-1) and j<(nlon_lr-1)):
				#print str(i).zfill(4),",",str(j).zfill(4)
				hsmth[i,j]=htmp[i,j]+smcf*(0.5*(htmp[i+1,j]+htmp[i-1,j])-htmp[i,j])
			else:
				hsmth[i,j]=ds.HGT_M[i,j]	

	ds_smth=xa.Dataset(
				data_vars={
					'HGT_M':(['south_north','west_east'],hsmth),
					}, 
					coords={'south_north':ds.south_north, 'west_east':ds.west_east}
				)	
	ds=ds.rename({'HGT_M':'HGT_O'})
	ds_new = ds.merge(ds_smth)
	
	print "  max. difference: {:4.1f} m".format(float((ds_new.HGT_O-ds_new.HGT_M).max().values))
	print "  std-dev        : {:4.1f} m".format(float((ds_new.HGT_O-ds_new.HGT_M).std().values))
	
	return ds_new


# ----------------------------------------------------------------------
# load_hr_dem(source, script_path)
#
# short: loads the high resolution DEM depending on the resolution
#		 of the topography that is to be generated
#
# arguments:		
#	source	... which DEM to choose as source for upscaling
#	script_path ... path of act.py to allow the script to find
#					the data files.
#
# returns:
#	dahr	... the high resolution topography
#
# description:
#	depending on the resolution requested for the topography that
#	is top be generated by act.py, the high resolution DEM is chosen.
#	e.g. if the minimum horizontal resolution is greater then 1' then
#	ETOPO1_Ice may be used, otherwise a higher resolution DEM is
#	necessary (such as SRTM for instance)
#	This routine requires the corresponding DEMs to be accessible in the
#	paths currently defined in this routine by the variable
#	script_path.
#
# ----------------------------------------------------------------------
def load_hr_dem(source,script_path):
	if source=='etopo1':
		data_path = './data/'
		dahr      = xa.open_dataset(script_path+'/'+data_path+"/"+ETOPO1I_FN).Band1
	elif source=='srtm':                # currently only works for New Zealand. Could use pySRTM
		data_path = '.data/srtm/'
		hr0       = xa.open_dataset(script_path+'/'+data_path+'srtm_70_21.nc').Band1.fillna(0)
		hr1       = xa.open_dataset(script_path+'/'+data_path+'srtm_70_22.nc').Band1.fillna(0)
		hr2       = xa.open_dataset(script_path+'/'+data_path+'srtm_71_21.nc').Band1.fillna(0)
		hr3       = xa.open_dataset(script_path+'/'+data_path+'srtm_71_22.nc').Band1.fillna(0)

		dahr      = hr0.combine_first(hr1).combine_first(hr2).combine_first(hr3)
	else:
		print " unknown data source: {:s}".format(source)
		
	return dahr

# set all potential options to initial values
lon1	= None
lon2    = None
lonc   	= None
lat0  	= None
lat1  	= None
lat2  	= None
dx    	= None
dy    	= None
nx		= None
ny		= None
name  	= None
smooth	= False
preview = False
plots	= False

try:
	opts, args = getopt.getopt(sys.argv[1:],"",["lonc=","lat0=","lat1=","lat2=","dx=","dy=","lon1=","lon2=","nx=","ny=","name=","smooth","preview","plots"])
except getopt.GetoptError, exc:
	print exc.msg
	print 'error occured, problem with parameters!'
	sys.exit(2)
for opt, arg in opts:
	if opt in ('--lonc'):
		lonc	= float(arg)
	elif opt in ('--lon1'):
		lon1   = float(arg)
	elif opt in ('--lon2'):
		lon2   = float(arg)
	elif opt in ('--lat0'):
		lat0 	= float(arg)
	elif opt in ('--lat1'):
		lat1	= float(arg)
	elif opt in ('--lat2'):
		lat2	= float(arg)
	elif opt in ('--dx'):
		dx   	= float(arg)
	elif opt in ('--dy'):
		dy   	= float(arg)
	elif opt in ('--nx'):
		nx   	= int(arg)
	elif opt in ('--ny'):
		ny   	= int(arg)
	elif opt in ('--name'):
		name 	= arg
	elif opt in ('--smooth'):
		smooth  = True
	elif opt in ('--preview'):
		preview = True
	elif opt in ('--plots'):
		plots   = True
		

none_pars_str=""
if (lonc is None) and (lon1 is None) and (lon2 is None):
	none_pars_str+="lonc or lon1 and lon2 "
if lon1 is None and not (lon2 is None):
	none_pars_str+='lon1 '
if lon2 is None and not (lon1 is None):
	none_pars_str+='lon2'
if (lat0 is None):
	none_pars_str=none_pars_str+"lat0 "
if (lat1 is None):
	none_pars_str=none_pars_str+"lat1 "
if (lat2 is None):
	none_pars_str=none_pars_str+"lat2 "
if (dy is None):
	none_pars_str+='dy '
if (dx is None):
	none_pars_str+='dx '

if (name is None):
	none_pars_str=none_pars_str+"name "
if (none_pars_str!=""):
	print " missing the following parameters: {:s}".format(none_pars_str)
	sys.exit(1)

# try to estimate nx and ny
if (nx is None) and (ny is None):
	dlon_km  = haversine(lon1,np.min([lat1,lat2]),lon2,np.min([lat1,lat2]))
	nx_guess = int(dlon_km*1000.0/dx)
	dlat_km  = haversine(lon1,lat1,lon1,lat2)
	ny_guess = int(dlat_km*1000.0/dy)
	nx       = 1.7*nx_guess
	ny       = 1.7*ny_guess						
	lonc	 = 0.5*(lon1+lon2)
	print '    estimated nx and ny as {:n} and {:n}'.format(nx,ny)
else:
	print '    set nx and ny to {:n} and {:n}'.format(nx,ny)
#sys.exit(1)


# create the correct projection
x0      = -4e5
y0      = -5e5
pwrf	= '+proj=lcc +lat_1={:2.5f} +lat_2={:2.5f} ' \
		  '+lat_0={:2.5f} +lon_0={:2.5f} ' \
          '+x_0={:f} +y_0={:f}'.format(lat1,lat2,lat0,lonc,0,0)
grid	= salem.Grid(nxny=(nx, ny), x0y0=(x0, y0), dxdy=(dx, dy), proj=pwrf)

lon1d	= np.min(grid.ll_coordinates[0])
lon2d	= np.max(grid.ll_coordinates[0])
lat1d	= np.min(grid.ll_coordinates[1])
lat2d	= np.max(grid.ll_coordinates[1])

print '  * {:3.1f}N {:3.1f}W to {:3.1f}N {:3.1f}W at {:2.1f}m x {:2.1f}m as {:s}'.format(lon1d,lat1d,lon2d,lat2d,dx,dy,name)

# output an overview of the topography
sm     = salem.Map(grid)
if preview:
	sm.visualize();
	print '    creating preview file...'
	plt.savefig('./{:s}_preview.pdf'.format(name))
	sys.exit(0)
	
script_path = os.path.dirname(os.path.realpath(__file__))


# calculate an average dlon and dlat to decided which data source to use
# equation from arclength at latitude theta where r=R*cos(theta)
# and b = phi * r
dlon0 = ((dx/1000.0)/(R*np.cos(d2r(lat0))))*(180.0/np.pi)*60.0
dlon1 = ((dx/1000.0)/(R*np.cos(d2r(lat1))))*(180.0/np.pi)*60.0
dlat  = ((dy/1000.0)/R)*(180.0/np.pi)*60.0	# from arclength b = phi * R than converted to degree and from there to minutes
dlon  = 0.5*(dlon0+dlon1)
dmin  = np.min([dlon,dlat])

# fall back on srtm dataset for higher resolutions
if (dmin < 1) or (dmin < 1):
	source = 'srtm'
else:
	source      = 'etopo1'
print '    minimum spacing is {:2.2f}'.format(dmin)
print '     => source dataset is {:s}...'.format(source)
topo_in     = load_hr_dem(source,script_path)
topo_in		= topo_in.sel(lon=slice(lon1d,lon2d),lat=slice(lat1d,lat2d))				# select the subset we need so that not the entire dataset is used
gr_hr 		= topo_in.salem.grid.regrid(factor=5)   									# define new grid with 5 times nx and 5 times ny gridpoints (
topo_hr		= gr_hr.map_gridded_data(topo_in, grid=topo_in.salem.grid,interp='nearest')	# interpolate to high resolution grid
count   	= np.where(topo_hr >= 0, 1, 0)												# count function that counts every datapoints above sea level - used
																						# to determine over how many gridpoints we average our lower resolution
																						# topography lateron

# generating lookup table and plot indicating over how many gridpoints averaging was performed
count_lr, lut 			= grid.lookup_transform(count, grid=gr_hr, method=np.sum, return_lut=True)
count_lr 				= count_lr * 1.
count_lr[count_lr == 0] = np.NaN

# output gridpoints count for every dx.dy cell on the map
if plots:
	f,ax = plt.subplots(1,1,dpi=200)
	sm.set_data(count_lr, grid)
	sm.visualize(ax=ax);
	plt.savefig('./{:s}_gridpoint_numbers.pdf'.format(name))

print '    upscaling...'
topo_lr					= grid.lookup_transform(topo_hr, grid=gr_hr, method=np.mean, lut=lut)
topo_lr[ topo_lr < 0.0]	= 0																		# everything below sea level is set to zero elevation

if plots:
	# output plot of dx.dy topography
	f,ax = plt.subplots(1,1,dpi=400)
	sm.set_data(topo_lr, grid)
	sm.visualize(ax=ax);
	plt.savefig('./{:s}_topography_original.pdf'.format(name))

# convert the data to the xarray format
topo_data = topo_lr.data

topo_data[np.isnan(topo_data)] = 0              # if this is not done, there are np.nan values in the nc file
lm_data                        = topo_data*0.0
lm_data[topo_data>0]           = 1.0            # lowest elevation is the sea level

topo_ds = xa.Dataset(data_vars={ 
            'HGT_M':(('south_north','west_east'),topo_data),
            'LANDMASK':(('south_north','west_east'),lm_data)
        },
                coords={
                        'XLONG_M': (('south_north','west_east'), grid.ll_coordinates[0]),
                        'XLAT_M': (('south_north','west_east'),  grid.ll_coordinates[1])
                       }
                )

if smooth == True:
	print '    filtering...'
	topo_ds=filter_topo(topo_ds)
print '    saving...'
topo_ds.to_netcdf('./'+name+'.nc')
print '  done'

sys.exit(0)

