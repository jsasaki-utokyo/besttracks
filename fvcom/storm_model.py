import os
import pickle
import PyFVCOM as pf
import netCDF4
from cftime import num2date, date2num
import numpy as np
from pyproj import Proj
from besttracks import parseJMA
import pandas as pd
import datetime
from dateutil.parser import parse

def read_fvcom_grd(fvcom_grid, utm='54'):
    '''
    Reading FVCOM grid
    
    Parameters
    ----------
    fvcom_grid : pickle or text
        FVCOM grid in UTM
    utm: int
        UTM zone number
    
    Returns
    -------
    x, y, xc, yc, nv, node, nele
    
    '''
    fvcom_grid_pickle = fvcom_grid[:-3] + "pickle"
    if os.path.isfile(fvcom_grid_pickle):
        ## If serialized file exists, it is to be opened.
        with open(fvcom_grid_pickle, "rb") as file:
            mesh = pickle.load(file)
            print("Reading the FVCOM grid from pickle. Ensure the serialized grid data is up-to-date.")
    else:
        print("Reading the FVCOM grid from", fvcom_grid, "and dumping into", fvcom_grid_pickle, 
              "for efficiency." )
        mesh = pf.grid.Domain(fvcom_grid, 'cartesian', utm)
        with open(fvcom_grid_pickle, "wb") as file:
            pickle.dump(mesh, file)
    x = mesh.grid.x
    y = mesh.grid.y
    xc = mesh.grid.xc
    yc = mesh.grid.yc
    nv = mesh.grid.triangles
    node = mesh.dims.node
    nele = mesh.dims.nele
    return x, y, xc, yc, nv, node, nele

def utm2geographic(x, y, xc, yc, zone=54, hemi='N'):
    '''
    Converting from UTM to Geographic coords
    
    Parameters
    ----------
    x, y, xc, yc
    zone: int
        UTM zone number
    hemi: 'N' or 'S'
        Nothern or Southern hemisphere
    
    Returns
    -------
    lon, lat, lonc, latc
    
    '''
    if hemi == 'S':
        y  = y  - 10000000
        yc = yc - 10000000
    ## Define coordinate converter p
    p = Proj(proj='utm', zone=zone, ellps='WGS84')
    lon, lat = p(x, y, inverse=True)
    lonc, latc = p(xc, yc, inverse=True)
    return lon, lat, lonc, latc 

def read_fvcom_nc(fvcom_nc):
    print("Reading the FVCOM grid from FVCOM output netCDF")
    FVCOM = netCDF4.Dataset(fvcom_nc, 'r')
    ## FVCOM = ncread(fvcom_nc)  ## too slow
    x=FVCOM.variables['x'][:]   ## nodes
    y=FVCOM.variables['y'][:]
    xc=FVCOM.variables['xc'][:] ## elements
    yc=FVCOM.variables['yc'][:]
    nv = FVCOM.variables['nv'][:].T - 1  ## tri elemnents; offset for Python indexing.
    node = np.arange(1, x.size+1)
    nele = np.arange(1, xc.size+1)
    return x, y, xc, yc, nv, node, nele

def besttrack(fpath, ID='201919', zone=54, sampling='H'):
    df = parseJMA(fpath)
    df = df[df.ID == ID]
    df = df.set_index('TIME')
    # df = df.tz_localize('UTC')  ## Set time zone as UTC
    df = df.resample(sampling).interpolate()  ## Up-sampling in every hour with interpolation
    p = Proj(proj='utm', zone=zone, ellps='WGS84')
    XUTM, YUTM = p(df.LON.values, df.LAT.values, inverse=False)  ## convert to UTM using Proj()
    df['XUTM'] = XUTM
    df['YUTM'] = YUTM
    ## Calculate translation speed
    ## The last row is a copy of the previous row so that the number of rows does not change. 
    TRANSLATION = np.sqrt(np.square(XUTM[1:]-XUTM[:-1]) + np.square(YUTM[1:]-YUTM[:-1])) / 3600
    TRANSLATION = np.append(TRANSLATION, TRANSLATION[-1])  ## Append the last row
    df['TRANSLATION'] = TRANSLATION
    df.fillna(method='ffill', inplace=True)
    return df

def days_since_calendar_date(timestamp, calendar_date='1858-11-17 00:00:00'):
    '''
    Calculate days in real since calendar_date in UTM.
    
    Parameters
    ----------
    timestamp: pandas Timestamp
    calendar_date: str
        Default modified julian day in FVCOM
    
    Returns
    -------
    days: real
    '''
    calendar_date = parse(calendar_date)
    #print('calendar_date = ', calendar_date)
    days_int = (timestamp - calendar_date).days
    within_day = ((timestamp - calendar_date).seconds / (24*3600))
    days = days_int + within_day
    return days

def netcdf_storm(ncfile, x, y, lon, lat, xc, yc, lonc, latc, nv, format="NETCDF4",
                 title = "FVCOM storm forcing file",
                 institution = "The University of Tokyo",
                 source = "FVCOM grid (unstructured) surface forcing",
                 references = "http://fvcom.smast.umassd.edu, http://codfish.smast.umassd.edu",
                 Conventions = "CF-1.0",
                 units = "days since 1858-11-17 00:00:00",
                 CoordinateSystem = "Cartesian",
                 CoordinateProjection = "init=epsg:32654"):
    '''
    Create and define netCDF for storm model contining wind velocity and surface pressure.
    
    Parameters
    ----------
    ncfile(str): netCDF file name
    x(ndarray): 1D real FVCOM node
    y(ndarray): 1D real FVCOM node
    lon(ndarray): 1D real FVCOM node
    lat(ndarray): 1D real FVCOM node
    xc(ndarray): 1D real FVCOM nele
    yc(ndarray): 1D real FVCOM nele
    lonc(ndarray): 1D real FVCOM nele (cell center)
    latc(ndarray): 1D real FVCOM nele (cell center)
    nv(ndarray): 2D int FVCOM nv(three, nele)
    global_atts(dict): dictionary for global attributes
    format(str): netCDF format: NETCDF3_64BIT_OFFSET(default), NETCDF4_CLASSIC, NETCDF4
    
    Returns
    -------
    netCDF4.Dataset
    '''

    nc = netCDF4.Dataset(ncfile, 'w', format=format)  ## Sometimes does not work with large data.
    # nc.set_fill_off()
    time = nc.createDimension("time", None)
    node = nc.createDimension("node", len(lat))
    nele  = nc.createDimension("nele", len(latc))
    three = nc.createDimension("three", 3)
    ## Coordinate variables
    times = nc.createVariable("time", "f8", ("time",))  ## f8 required; f4 causing inaccurate time
    nodes = nc.createVariable("node", "i4", ("node",))
    neles = nc.createVariable("nele", "i4", ("nele",))
    ## Variables
    xs = nc.createVariable("x", "f4", ("node",))
    ys = nc.createVariable("y", "f4", ("node",))
    lons = nc.createVariable("lon", "f4", ("node",))
    lats = nc.createVariable("lat", "f4", ("node",))
    xsc = nc.createVariable("xc", "f4", ("nele",))
    ysc = nc.createVariable("yc", "f4", ("nele",))
    lonsc = nc.createVariable("lonc", "f4", ("nele",))
    latsc = nc.createVariable("latc", "f4", ("nele",))
    nvs = nc.createVariable("nv", "i4", ("nele","three",))
    sp = nc.createVariable("air_pressure", "f4", ("time", "node",))
    u_w_x  = nc.createVariable("u_w_x",  "f4", ("time", "nele",))
    u_w_y  = nc.createVariable("u_w_y",  "f4", ("time", "nele",))
    ## Global attributes
    nc.title = title
    nc.institution = institution
    nc.source = source
    nc.history = "created " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    nc.references = references
    nc.Conventions = Conventions
    nc.units = units
    nc.CoordinateSystem = CoordinateSystem
    nc.CoordinateProjection = CoordinateProjection
    ## Coordinate variable attributes
    times.units = "days since 1858-11-17 00:00:00" ## Modified Julian Day
    times.calendar = "standard"
    ## Variable attributes
    xs.long_name = "nodal x-coordinate"
    xs.units = "meters"
    ys.long_name = "nodal y-coordinate"
    ys.units = "meters"
    lons.long_name = "nodal longitude"
    lons.standard_name = "longitude"
    lons.units = "degrees_east"
    lats.long_name = "nodal latitude"
    lats.standard_name = "latitude"
    lats.units = "degrees_north"
    xsc.long_name = "zonal x-coordinate"
    xsc.units = "meters"
    ysc.long_name = "zonal y-coordinae"
    ysc.units = "meters"
    lonsc.long_name = "zonal longitude"
    lonsc.standard_name = "longitude"
    lonsc.units = "degrees_east"
    latsc.long_name = "zonal latitude"
    latsc.standard_name = "latitude"
    latsc.units = "degrees_north"
    nvs.long_name = "nodes surrounding element"
    sp.units = "hPa"
    sp.long_name = "surface air pressure"
    u_w_x.units = "m/s"
    u_w_x.long_name = "x component of wind velocity"
    u_w_y.units = "m/s"
    u_w_y.long_name = "y component of wind velocity"
    ## Writing coordinate variable data
    #print(lat)
    xs[:] = x
    ys[:] = y
    lats[:] = lat
    lons[:] = lon
    xsc[:] = xc
    ysc[:] = yc
    latsc[:] = latc
    lonsc[:] = lonc
    nvs[:,:] = nv
    return nc

def parametric_model(df, ncfile, x, y, lon, lat, xc, yc, lonc, latc, nv,
                     P_0=1013.25, alpha=30, **kwargs):
    xt_c = df.XUTM.values  ## typhoon track center x in UTM
    yt_c = df.YUTM.values  ## typhoon track center y in UTM
    sp_c  = df.PRS.values   ## typhoon track central surface pressure in hPa
    translation_speed = df.TRANSLATION.values
    num_tracks = xt_c.size - 1 ## should be -1 to obtain typhoon translation vector (one more track required)
    timestamps = df.index
    ## Create netCDF
    nc = netcdf_storm(ncfile, x, y, lon, lat, xc, yc, lonc, latc, nv)
    u_w1_x = nc.createVariable("u_w1_x", "f4", ("time", "nele",))
    u_w1_y = nc.createVariable("u_w1_y", "f4", ("time", "nele",))
    u_w2_x = nc.createVariable("u_w2_x", "f4", ("time", "nele",))
    u_w2_y = nc.createVariable("u_w2_y", "f4", ("time", "nele",))
    u_w1_x.units = "m/s"
    u_w1_x.long_name = "x component of cyclonic wind velocity"
    u_w1_y.units = "m/s"
    u_w1_y.long_name = "y component of cyclonic wind velocity"
    u_w2_x.units = "m/s"
    u_w2_x.long_name = "x component of translation velocity"
    u_w2_y.units = "m/s"
    u_w2_y.long_name = "y component of translation velocity"

    ## Surface pressure of Myers model
    P_0 = P_0
    r_min = 10
    r_max = (1.633 * sp_c - 1471.35) * 1000  ## when sp_c >= 950 hPa
    r_max_tmp = (0.769 * sp_c - 650.55) * 1000  ## when 880 hPa <= sp_c < 950 hPa
    index_sp_c_lt_950, = np.where(sp_c < 950)  ## array index where sp_c < 950 hPa 
    r_max[index_sp_c_lt_950] = r_max_tmp[index_sp_c_lt_950]
    index_sp_c_lt_880, = np.where(sp_c < 880)
    #index_sp_c_lt_880.size
    if index_sp_c_lt_880.size > 0:
        formatted = f"WARNING: sp_c at indices of {index_sp_c_lt_880} < 880 hPa"
        print(formatted)

    ## sp: atm pressure at each node at i_track
    for i_track in np.arange(num_tracks):  ## loop for tracks
        r = np.sqrt(np.square(x - xt_c[i_track].item()) + np.square(y - yt_c[i_track].item()))
        sp = np.where(r > r_min, (sp_c[i_track]+(P_0-sp_c[i_track]) * np.exp(-r_max[i_track]/r)),
                      sp_c[i_track])
        nc.variables["air_pressure"][i_track, :] = sp
        nc.variables["time"][i_track] = days_since_calendar_date(timestamps[i_track])

    ## Wind field of Mitsuta-Fujii's model (1987)
    ## Referred to [Tajima et al. (2016)](https://doi.org/10.1142/S0578563416400027)
    ## `x_c`, `y_c` (`latc`): at elements.
    ## `xt_c` and `yt_c`: at typhoon track center
    C1=2/3
    C2=2/3
    rho_a=1.225
    omega = 2 * np.pi/(24*3600) ## angular frequency of earth rotation
    alpha = alpha ## anticlockwise angle in dgree from the tangential direction of wind

    alpha = np.pi*alpha/180  ## degree -> radian
    for i_track in np.arange(num_tracks):  ## loop for typhoon tracks
        r = np.sqrt(np.square(xc - xt_c[i_track]) + np.square(yc - yt_c[i_track]))
        r = np.where(r >= r_min, r, r_min)  ## check: originally np.where(r >= r_min, r-r_min, r_min) but r-r_min may be zero. 
        f_c = 2 * omega * np.sin(np.pi*latc/180)
        dpa_dr = (P_0 - sp_c[i_track]) * 100 * r_max[i_track]/np.square(r) * np.exp(-r_max[i_track]/r)
        u_w1 = C1 * (-f_c*r/2 + np.sqrt(np.square(f_c*r/2) + r/rho_a * dpa_dr))
        u_gr = C1 * (-f_c*r/2 + np.sqrt(np.square(f_c*r/2)
             + 100 * (P_0 - sp_c[i_track])*r_max[i_track]/(rho_a*r) * np.exp(-r_max[i_track]/r) ))
        dpa_dr_r_max = (P_0 - sp_c[i_track]) * 100 * np.exp(-1)/r_max[i_track]  ## multply 100 to convert from hPa to Pa
        ## tangential direction
        u_w1_x_tmp = -(yc - yt_c[i_track])/r * u_w1  # x-component
        u_w1_y_tmp =  (xc - xt_c[i_track])/r * u_w1  # y-component
        ## rotating 30 degrees anticlockwise from the tangential direction
        u_w1_x = u_w1_x_tmp * np.cos(alpha) - u_w1_y_tmp * np.sin(alpha)
        u_w1_y = u_w1_x_tmp * np.sin(alpha) + u_w1_y_tmp * np.cos(alpha)
    
        u_w1_r_max = C1*(-f_c*r_max[i_track]/2 + np.sqrt(np.square(f_c*r_max[i_track]/2)
                   + r_max[i_track]/rho_a * dpa_dr_r_max))
        u_w2 = C2 * u_w1/u_w1_r_max * translation_speed[i_track]
    
        dr_t = np.sqrt(np.square(xt_c[i_track+1]-xt_c[i_track]) + np.square(yt_c[i_track+1]-yt_c[i_track]))
        u_w2_x = (xt_c[i_track+1] - xt_c[i_track])/dr_t * u_w2
        u_w2_y = (yt_c[i_track+1] - yt_c[i_track])/dr_t * u_w2
    
        u_w_x = u_w1_x + u_w2_x
        u_w_y = u_w1_y + u_w2_y

        nc.variables["u_w1_x"][i_track, :] = u_w1_x
        nc.variables["u_w1_y"][i_track, :] = u_w1_y
        nc.variables["u_w2_x"][i_track, :] = u_w2_x
        nc.variables["u_w2_y"][i_track, :] = u_w2_y
        nc.variables["u_w_x"][i_track, :] = u_w_x
        nc.variables["u_w_y"][i_track, :] = u_w_y
    nc.close()

def era_model(era_uv_nc, era_sp_nc, ncf_era, x, y, lon, lat, xc, yc, lonc, latc, nv):
    ## The following results in an error because of the conflict of nc.
    ## nc = netcdf_storm(ncf_era, x, y, lon, lat, xc, yc, lonc, latc, nv)
    with netCDF4.Dataset(era_uv_nc, 'r') as nc:
        lon_era_uv = nc.variables['longitude'][:]
        lat_era_uv = nc.variables['latitude'][:]
        time_era_uv = nc.variables['time'][:]
        u_era = nc.variables['u10'][:]
        v_era = nc.variables['v10'][:]
    num_lon_era_uv = lon_era_uv.size
    num_lat_era_uv = lat_era_uv.size
    num_time_era_uv = time_era_uv.size

    with netCDF4.Dataset(era_sp_nc, 'r') as nc:
        lon_era_sp = nc.variables['longitude'][:]
        lat_era_sp = nc.variables['latitude'][:]
        time_era_sp = nc.variables['time'][:]
        sp_era = nc.variables['sp'][:]
        ## Check whether both times are identical, or raise an error.
        if np.all(time_era_uv == time_era_sp):
            print("time_era_uv == time_era_sp; use time_era")
            time_era = time_era_uv
            time_era_units = nc.variables['time'].units
            time_era_calendar = nc.variables['time'].calendar
            print(time_era_units)
            print(time_era_calendar)
        else:
            raise ValueError("ERROR: time_era /= time_era_uv")
    num_lon_era_sp = lon_era_sp.size
    num_lat_era_sp = lat_era_sp.size
    num_time_era_sp = time_era_sp.size
    
    ## Find a ERA grid point `(i,j)` where the grid area between `(i,j)` and `(i+1,j+1)`
    ## contains an FVCOM node or cell center
    '''
        i       i+1    *: FVCOM grid node or cell center index n
    j   |--------|     i: Longitude grid index of ERA   
        |  * n   |     j: Latitude grid indes of ERA (positive downward)
    j+1 |--------|
    '''
    ## Interpolate ERA surface pressure at FVCOM node
    dlon_era_sp = lon_era_sp[1] - lon_era_sp[0]
    dlat_era_sp = lat_era_sp[0] - lat_era_sp[1]
    i_sp = (np.trunc((lon - lon_era_sp[0]) / dlon_era_sp)).astype(int)
    j_sp = (np.trunc((lat_era_sp[0] - lat) / dlat_era_sp)).astype(int)
    cx0 = (lon - lon_era_sp[i_sp]) / (lon_era_sp[i_sp+1] - lon_era_sp[i_sp])
    cx1 = (lon_era_sp[i_sp + 1] - lon) / (lon_era_sp[i_sp + 1] - lon_era_sp[i_sp])
    cy0 = (lat_era_sp[j_sp] - lat) / (lat_era_sp[j_sp] - lat_era_sp[j_sp + 1])
    cy1 = (lat - lat_era_sp[j_sp + 1]) / (lat_era_sp[j_sp] - lat_era_sp[j_sp + 1])
    dates_era = num2date(time_era, units=time_era_units, calendar=time_era_calendar)
    days_era = date2num(dates_era, units="days since 1858-11-17 00:00:00", calendar="gregorian")
    ## ERAの時間とFVCOMの時間を合わせ，netCDFはFVCOMの時間で作る
    ## Matching ERA time and FVCOM time; creating netCDF time using FVCOM time
    nc = netcdf_storm(ncf_era, x, y, lon, lat, xc, yc, lonc, latc, nv)
    for time in np.arange(time_era.size):
        sp_era_fvcom_j0 = cx1 * sp_era[time, j_sp, i_sp] + cx0 * sp_era[time, j_sp, i_sp+1]
        #print(sp_era_fvcom_j0.shape)
        sp_era_fvcom_j1 = cx1 * sp_era[time, j_sp+1, i_sp] + cx0 * sp_era[time, j_sp+1, i_sp+1]
        #print(sp_era_fvcom_j1.shape)
        sp_era_fvcom = cy1 * sp_era_fvcom_j0 + cy0 * sp_era_fvcom_j1
        #print(sp_era_fvcom.shape)
        nc.variables["time"][time] = days_era[time]
        nc.variables["air_pressure"][time, :] = sp_era_fvcom

    ## Interpolate ERA wind field at FVCOM cell center
    dlon_era_uv = lon_era_uv[1] - lon_era_uv[0]
    dlat_era_uv = lat_era_uv[0] - lat_era_uv[1]
    i_uv = (np.trunc((lonc - lon_era_uv[0]) / dlon_era_uv)).astype(int)
    j_uv = (np.trunc((lat_era_uv[0] - latc) / dlat_era_uv)).astype(int)
    cx0 = (lonc - lon_era_uv[i_uv]) / (lon_era_uv[i_uv+1] - lon_era_uv[i_uv])
    cx1 = (lon_era_uv[i_uv + 1] - lonc) / (lon_era_uv[i_uv + 1] - lon_era_uv[i_uv])
    cy0 = (lat_era_uv[j_uv] - latc) / (lat_era_uv[j_uv] - lat_era_uv[j_uv + 1])
    cy1 = (latc - lat_era_uv[j_uv + 1]) / (lat_era_uv[j_uv] - lat_era_uv[j_uv + 1])
                                                                      
    for time in np.arange(time_era.size):
        u_era_fvcom_j0 = cx1 * u_era[time, j_uv, i_uv] + cx0 * u_era[time, j_uv, i_uv+1]
        #print(sp_era_fvcom_j0.shape)
        u_era_fvcom_j1 = cx1 * u_era[time, j_uv+1, i_uv] + cx0 * u_era[time, j_uv+1, i_uv+1]
        #print(sp_era_fvcom_j1.shape)
        u_era_fvcom = cy1 * u_era_fvcom_j0 + cy0 * u_era_fvcom_j1
        #print(sp_era_fvcom.shape)
        nc.variables["u_w_x"][time, :] = u_era_fvcom

        v_era_fvcom_j0 = cx1 * v_era[time, j_uv, i_uv] + cx0 * v_era[time, j_uv, i_uv+1]
        v_era_fvcom_j1 = cx1 * v_era[time, j_uv+1, i_uv] + cx0 * v_era[time, j_uv+1, i_uv+1]
        v_era_fvcom = cy1 * v_era_fvcom_j0 + cy0 * v_era_fvcom_j1
        nc.variables["u_w_y"][time, :] = v_era_fvcom
    nc.close()

def hybrid_model(df, ncf_pm, ncf_era, ncf_hv, x, y, lon, lat, xc, yc, lonc, latc, nv):
    ''' Open parametric model netCDF and ERA netCDF
        Extract the start (`ids`) and end (`ide`) indices of the common time range for
        PRM (parametric) and ERA
        Extract the same time range (from `datetime_s` to `datetime_e`) for besttrack DataFrame
    '''
    nc = netcdf_storm(ncf_hv, x, y, lon, lat, xc, yc, lonc, latc, nv)

    PRM = netCDF4.Dataset(ncf_pm, "r")
    ERA = netCDF4.Dataset(ncf_era, "r")
    time_PRM = PRM.variables["time"]
    time_ERA = ERA.variables["time"]
    # print(num2date(time_PRM[:], units=time_PRM.units, calendar=time_PRM.calendar))
    # print(num2date(time_ERA[:], units=time_ERA.units, calendar=time_ERA.calendar))
    day_s, day_e = max(time_PRM[0], time_ERA[0]).item(), min(time_PRM[-1], time_ERA[-1]).item()
    print(day_s, day_e)
    ids_PRM, ide_PRM = np.where(time_PRM[:] == day_s)[0].item(), np.where(time_PRM[:] == day_e)[0].item()
    print(ids_PRM, ide_PRM)
    ids_ERA, ide_ERA = np.where(time_ERA[:] == day_s)[0].item(), np.where(time_ERA[:] == day_e)[0].item()
    print(ids_ERA, ide_ERA)
    datetime_s = str(num2date(time_PRM[ids_PRM], units=time_PRM.units, calendar=time_PRM.calendar))
    datetime_e = str(num2date(time_PRM[ide_PRM], units=time_PRM.units, calendar=time_PRM.calendar))

    ## Set inner (`r_i`) and outer (`r_o`) radii and band width (`b_w = r_o - r_i`) in meters
    r_i = 300000; r_o = 500000; b_w = r_o - r_i
    n = -1
    for n_PRM in np.arange(ide_PRM-ids_PRM+1):
        n += 1
        ## Gets ERA time index (n_ERA) corresponding to 0-indexed PRM's n
        n_ERA = n_PRM + ids_ERA - ids_PRM
        datetime_BST = str(num2date(time_PRM[n_PRM], units=time_PRM.units, calendar=time_PRM.calendar))
        xt_c = df.loc[datetime_BST, 'XUTM']  ## x (UTM) of typhoon center
        yt_c = df.loc[datetime_BST, 'YUTM']  ## y (UTM) of typhoon center
        dist_sp = np.sqrt(np.square(x  - xt_c) + np.square(y  - yt_c))
        dist_uv = np.sqrt(np.square(xc - xt_c) + np.square(yc - yt_c))

        sp = ERA['air_pressure'][n_ERA, :]
        sp = np.where(dist_sp <= r_i, PRM['air_pressure'][n_PRM, :], sp)
        ## dist_sp = dist_sp when FVCOM node is in the region between inner and outer circles,
        #  otherwise = np.nan
        ## dist_sp between inner and outer circles
        dist_sp = np.where((dist_sp > r_i) & (dist_sp < r_o), dist_sp, np.nan)
        dr_i = dist_sp - r_i  ## distance between FVCOM True node and inner circle
        dr_o = r_o - dist_sp  ## distance between FVCOM True node and outer circle
        sp = np.where(~np.isnan(dist_sp), dr_i/b_w * ERA['air_pressure'][n_ERA,:]
                      + dr_o/b_w * PRM['air_pressure'][n_PRM,:], sp)

        nc.variables["time"][n] = PRM["time"][n_PRM]
        nc.variables["air_pressure"][n, :] = sp

        u_w_x = ERA['u_w_x'][n_ERA, :]
        u_w_y = ERA['u_w_y'][n_ERA, :]
        u_w_x = np.where(dist_uv <= r_i, PRM['u_w_x'][n_PRM, :], u_w_x)
        u_w_y = np.where(dist_uv <= r_i, PRM['u_w_y'][n_PRM, :], u_w_y)
        ## dist_uv = dist_uv when FVCOM nele is in the region between inner and outer circles,
        #  otherwise = np.nan
        ## dist_uv between inner and outer circles
        dist_uv = np.where((dist_uv > r_i) & (dist_uv < r_o), dist_uv, np.nan)
        dr_i = dist_uv - r_i  ## distance between FVCOM True nele and inner circle
        dr_o = r_o - dist_uv  ## distance between FVCOM True nele and outer circle
        u_w_x = np.where(~np.isnan(dist_uv), dr_i/b_w * ERA['u_w_x'][n_ERA,:] + dr_o/b_w * PRM['u_w_x'][n_PRM,:], u_w_x)
        u_w_y = np.where(~np.isnan(dist_uv), dr_i/b_w * ERA['u_w_y'][n_ERA,:] + dr_o/b_w * PRM['u_w_y'][n_PRM,:], u_w_y)
        nc.variables["u_w_x"][n, :] = u_w_x
        nc.variables["u_w_y"][n, :] = u_w_y
    PRM.close()
    ERA.close()
    nc.close()

def write_fvcom_forcing(ncfile, nele, node, global_attributes, nv, Times, x, y, lon, lat, 
                        xc, yc, lonc, latc, uwind_speed, vwind_speed, air_pressure,
                        format='NETCDF4'):
    '''
    
    '''
    dimensions = {'nele': nele, 'node': node, 'three': 3, 'time': None, 'DateStrLen': 26}
    with pf.preproc.WriteForcing(filename=ncfile, dimensions=dimensions,
                                 global_attributes=global_attributes, format=format) as nc:
        nc.write_fvcom_time(time=Times)
        nc.add_variable('x', x, ('node',), {'long_name': 'nodal x-coordinate', 'units': 'meters'})
        nc.add_variable('y', y, ('node',), {'long_name': 'nodal y-coordinate', 'units': 'meters'})
        nc.add_variable('lon', lon, ('node',), {'long_name': 'nodal longitude',
                        'standard_name': 'longitude', 'units': 'degrees_east'})
        nc.add_variable('lat', lat, ('node',), {'long_name': 'nodal latitude',
                        'standard_name': 'latitude', 'units': 'degrees_north'})
        nc.add_variable('xc', xc, ('nele',), {'long_name': 'zonal x-coordinate', 'units': 'meters'})
        nc.add_variable('yc', yc, ('nele',), {'long_name': 'zonal y-coordinate', 'units': 'meters'})
        nc.add_variable('lonc', lonc, ('nele',), {'long_name': 'zonal longitude', 
                        'standard_name': 'longitude', 'units': 'degrees_east'})
        nc.add_variable('latc', latc, ('nele',), {'long_name': 'zonal latitude',
                        'standard_name': 'latitude', 'units': 'degrees_north'})
        nc.add_variable('nv', nv.T, ('three', 'nele',),
                        {'long_name': 'nodes surrounding element'}, 'i')
        nc.add_variable('uwind_speed', uwind_speed, ('time', 'nele',),
                        {'long_name': 'Eastward Wind Speed', 'standard_name': 'Wind Speed',
                         'units': 'm/s', 'grid': 'fvcom_grid', 'type': 'data'})
        nc.add_variable('vwind_speed', vwind_speed, ('time', 'nele',),
                        {'long_name': 'Northward Wind Speed', 'standard_name': 'Wind Speed',
                         'units': 'm/s', 'grid': 'fvcom_grid', 'type': 'data'})
        nc.add_variable('air_pressure', air_pressure, ('time', 'node',),
                        {'long_name': 'surface air pressure', 'units': 'Pa', 'grid': 'fvcom_grid',
                         'coordinates': 'FVCOM Cartesian coordinates', 'type': 'data'})

def create_storm_models(casename, type=None, f_besttrack=None, ID=None, era_uv_nc=None, era_sp_nc=None, 
                        zone=54, sampling='H', P_0=1013.25, alpha=30,
                        ncf_pm='parametric_model.nc', ncf_era='era_model.nc', ncf_hv='hybrid_model.nc'):
    fvcom_grid = casename + "_grd.dat"
    x, y, xc, yc, nv, node, nele = read_fvcom_grd(fvcom_grid, utm=str(zone))
    lon, lat, lonc, latc = utm2geographic(x, y, xc, yc, zone=str(zone))
    if type is None:
        print('Creating Parametric, ERA, and Hybrid models')
        df = besttrack(f_besttrack, ID, zone=zone, sampling=sampling)
        parametric_model(df, ncf_pm, x, y, lon, lat, xc, yc, lonc, latc, nv, P_0=P_0, alpha=alpha)
        era_model(era_uv_nc, era_sp_nc, ncf_era, x, y, lon, lat, xc, yc, lonc, latc, nv)
        hybrid_model(df, ncf_pm, ncf_era, ncf_hv, x, y, lon, lat, xc, yc, lonc, latc, nv)
    elif type[0] == 'P':
        print('Creating Parametric model only')
        df = besttrack(f_besttrack, ID, zone=zone, sampling=sampling)
        parametric_model(df, ncf_pm, x, y, lon, lat, xc, yc, lonc, latc, nv, P_0=P_0, alpha=alpha)
    elif type[0] == 'E':
        print('Creating ERA model only')
        era_model(era_uv_nc, era_sp_nc, ncf_era, x, y, lon, lat, xc, yc, lonc, latc, nv)
    else:
        print("ERROR: No such type. Type should be None, 'PARAMETRIC', or 'ERA']")

def read_storm_model(ncfile, item, geo=False, format='NETCDF4'):
    with netCDF4.Dataset(ncfile, "r", format=format) as nc:
        var = nc.variables[item][:]
        node = len(nc.dimensions['node'])
        nele = len(nc.dimensions['nele'])
        time = nc.variables['time']
        time = pd.to_datetime(num2date(time[:], units=time.units, calendar=time.calendar, 
                              only_use_cftime_datetimes=False))
        if len(var[0,:]) == node:
            if geo:
                x = nc.variables['lon'][:]
                y = nc.variables['lat'][:]
            else:
                x = nc.variables['x'][:]
                y = nc.variables['y'][:]
            return time, x, y, var
        else:
            if geo:
                xc = nc.variables['lonc'][:]
                yc = nc.variables['latc'][:]
            else:
                xc = nc.variables['xc'][:]
                yc = nc.variables['yc'][:]
            return time, xc, yc, var
