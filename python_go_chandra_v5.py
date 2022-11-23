# coding: utf-8

# The code takes the corrected file from *sso_freeze* (hardwired by user) and peforms a corrdinate transformation on the X-ray emission to wrap the PSF around Jupiter

# Authors: Dale Weigt (D.M.Weigt@soton.ac.uk), apadpted from Randy Gladstone's 'gochandra' IDL script. Other contributors to the current (and/or previous) iteration(s) of the code are: Hunter Waite, Kurt Franke, Peter Ford, Seán McEntee, Caitríona Jackman, Will Dunn, Brad Snios, Ron Elsner, Ralph Kraft, and Graziella Branduardi-Raymont 

"""All the relevant packages are imported for code below"""

import go_chandra_analysis_tools as gca_tools # import the defined functions to analysis Chandra data nad perfrom coordinate transformations

import numpy as np
import pandas as pd
import scipy
from scipy import interpolate
from astropy.io import ascii
from astropy.io import fits as pyfits
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import colors
import matplotlib.gridspec as gridspec
import os
from astropy.time import Time

# Reading in configuration file containing any hard wired inputs.
import configparser
config = configparser.ConfigParser()
config.read('config.ini')


# AU to meter conversion - useful later on (probably a function built in already)
AU_2_m = 1.49598E+11
AU_2_km = 1.49598E+8

obsID = str(config['inputs']['obsID'])

# Accounting for different filepaths of ObsIDs that originlly had SAMP values and others that did not.
# df = pd.read_csv('ObsIDs_with_samp.txt', header=None, delimiter='\t')
# samp_ids = np.array(df.iloc[:,0])

# Reading in Chandra Event file, extracting all the relevant info and defining assumptions used in analysis <br>

# User is prompted to enter the file path of the corrected event file. The script finds the file from the selected folder and reads in all the relevent headers. The asusmptions used for the mapping are also defined here.

#Read in fits file
# if int(obsID) in samp_ids:
#     folder_path = '/Users/mcentees/Desktop/Chandra/' + str(obsID) + '/primary'
# else:
#     folder_path = '/Users/mcentees/Desktop/Chandra/' + str(obsID) + '/repro'

folder_path = str(config['inputs']['folder_path'])
cor_evt_location = []

# Script then searches through the folder looking the filename corresponding to the corrected file
for file in os.listdir(str(folder_path)):
    if file.startswith("hrcf") and file.endswith("pytest_evt2.fits"):
        cor_evt_location.append(os.path.join(str(folder_path), file))

# File is then read in with relevant header information extracted:
hdulist = pyfits.open(cor_evt_location[0], dtype=float)
matplotlib.rcParams['agg.path.chunksize'] = 10000
img_events=hdulist['EVENTS'].data # the data of the event file
img_head = hdulist[1].header # the header information of the event file
bigtime = img_events['time'] # time
bigxarr = img_events['X'] # x position of photons
bigyarr = img_events['Y'] # y position of photons
bigchannel = img_events['pha'] # pha channel the photons were found in
obs_id = img_head['OBS_ID'] # observation id of the event
date_start = img_head['DATE-OBS'] # Start date of observation
date_end = img_head['DATE-END'] # End date of observation
tstart = img_head['TSTART'] # the start and...
tend = img_head['TSTOP'] #... end time of the observation
sumamps = img_events['sumamps'] # reading in sumamps figure
samp = img_events['samp'] # reading in samp figure
pi_cal = img_events['pi']

cxo_tstart = str(Time(date_start, format='isot').decimalyear).replace('.','D')
cxo_tend = str(Time(date_end, format='isot').decimalyear).replace('.','D')

# reading in amplifier signals 
av1 = img_events['av1']
av2 = img_events['av2']
av3 = img_events['av3']

au1 = img_events['au1']
au2 = img_events['au2']
au3 = img_events['au3']

amp_sf = img_events['amp_sf'] # reading in amplifier scaling factor

# The date of the observation is read in...
datestart = img_head['DATE-OBS']
evt_date = pd.to_datetime(datestart) #... and coverted to datetiem format to allow the relevant information to be read to...
evt_hour = evt_date.hour
evt_doy = evt_date.strftime('%j')
evt_mins = evt_date.minute
evt_secs = evt_date.second
evt_DOYFRAC = gca_tools.doy_frac(float(evt_doy), float(evt_hour), float(evt_mins), float(evt_secs)) #... calculated a fractional Day of 
# Year (DOY) of the observation

ra_centre, ra_centre_rad = img_head['RA_NOM'], np.deg2rad(img_head['RA_NOM']) # the RA of Jupiter at the centre of the chip is read in as...
dec_centre, dec_centre_rad = img_head['DEC_NOM'], np.deg2rad(img_head['DEC_NOM']) #... well as Jupitr's DEC
j_rotrate = np.rad2deg(1.758533641E-4) # Jupiter's rotation period

hdulist.close()

# Assumptions used for mapping:
scale = 0.13175 # scale used when observing Jupiter using Chandra - in units of arcsec/pixel
fwhm = 0.8 # FWHM of the HRC-I point spread function (PSF) - in units of arcsec
psfsize = 25 # size of PSF used - in units of arcsec
alt = 400 # altitude where X-ray emission is assumed to occur in Jupiter's ionosphere - in units of km

# Alogrithm uses the start and end date from the observation to generate an epheremis file (from the JPL Horizons server) to use for analysis. The ephermeris file used takes CXO as the observer

"""Brad's horizons code to extract the ephemeris file"""

from astropy.time import Time                   #convert between different time coordinates
from astropy.time import TimeDelta              #add/subtract time intervals 
from astroquery.jplhorizons import Horizons     #automatically download ephemeris 

# The start and end times are taken from the horizons file.
tstart_eph=Time(tstart, format='cxcsec')
tstop_eph=Time(tend, format='cxcsec')
dt = TimeDelta(0.125, format='jd')

# Below sets the parameters of what observer the ephemeris file is generated form. For example, '500' = centre of the Earth, '500@-151' = CXO
obj = Horizons(id=599,location='500@-151',epochs={'start':tstart_eph.iso, 'stop':(tstop_eph+dt).iso, 'step':'1m'}, id_type='majorbody')
eph_jup = obj.ephemerides()

# Extracts relevent information needed from ephermeris file
cml_spline_jup = scipy.interpolate.UnivariateSpline(eph_jup['datetime_jd'], eph_jup['PDObsLon'],k=1)
lt_jup = eph_jup['lighttime']
sub_obs_lon_jup = eph_jup['PDObsLon']
sub_obs_lat_jup = eph_jup['PDObsLat']

# Adding angular diameter from JPL Horizons to use later to define radius of circular region within which photons are kept
ang_diam = max(eph_jup['ang_width'])

# Also adding tilt angle of Jupiter with respect to true North Pole
tilt_ang = np.mean(eph_jup['NPole_ang'])

# saving angular diameter and tilt angle in text file in order to plot ellipse in post-processing
np.savetxt(str(folder_path) + f'/{obsID}_JPL_ellipse_vals.txt', np.c_[ang_diam, tilt_ang], delimiter=',', header='angular diameter (arcsec),tilt angle (deg)', fmt='%s')

eph_dates = pd.to_datetime(eph_jup['datetime_str'])
eph_dates = pd.DatetimeIndex(eph_dates)
eph_doy = np.array(eph_dates.strftime('%j')).astype(int)
eph_hours = eph_dates.hour
eph_minutes = eph_dates.minute
eph_seconds = eph_dates.second

eph_DOYFRAC_jup = gca_tools.doy_frac(eph_doy, eph_hours, eph_minutes, eph_seconds) # DOY fraction from ephermeris data

jup_time = (eph_DOYFRAC_jup - evt_DOYFRAC)*86400.0 + tstart # local tiem of Jupiter

# Select Region for analysis

# Plots the photons (x,y) position on a grid of defined size in arcseconds (defualted at [-50,50] in both x and y). Jupiter is centred on the HRC instrument. The photon information form the defined 

# converting the x and y coordinates from the event file into arcseconds
bigxarr_region = (bigxarr - 16384.5)*0.13175
bigyarr_region = (bigyarr - 16384.5)*0.13175

# storing all photon data in text file - need this to calculate area for samp distributions later on
# np.savetxt(str(folder_path) + r"/%s_all_photons.txt" % obs_id, np.c_[bigxarr_region, bigyarr_region, bigtime, bigchannel, samp, sumamps, pi_cal, amp_sf, av1, av2, av3, au1, au2, au3])


# define the x, y, and pha channel limits 
xlimits, ylimits = [-50,50], [-50,50]
cha_min = 0
cha_max = max(bigchannel)

# the photon data is stored in a pandas dataframe 
evt_df = pd.DataFrame({'time': bigtime, 'x': bigxarr, 'y': bigyarr, 'pha': bigchannel})

# defines the region the photons will be selected from
indx = gca_tools.select_region(xlimits[0], xlimits[1],ylimits[0], ylimits[1],bigxarr_region,bigyarr_region,bigchannel,cha_min,cha_max)

# find the x and y position of the photons
x_ph = bigxarr_region[indx]
y_ph = bigyarr_region[indx]

# plots the selected region (sanity check: Jupiter should be in the centre)
fig, axes=plt.subplots(figsize=(7,7))
axes = plt.gca()


plt.plot(x_ph,y_ph, 'o', markersize=0.5,linestyle='None',color='blue')
plt.title('Selected Region (ObsID %s)' % obs_id)
plt.xlim(xlimits)
plt.ylim(ylimits)
print('')
print('')
print('Once you are happy with the selected region, close the figure window to continue analysis')
print('')
print('')
plt.show()

# saves the selected region as a text file
np.savetxt(str(folder_path) + r"/%s_selected_region_ellipse.txt" % obs_id, np.c_[x_ph, y_ph, bigtime[indx], bigchannel[indx], samp[indx], sumamps[indx], pi_cal[indx], amp_sf[indx], av1[indx], av2[indx], av3[indx], au1[indx], au2[indx], au3[indx]])

ph_data = ascii.read(str(folder_path) + r"/%s_selected_region_ellipse.txt" % obs_id) # read in the selected region data and...
ph_time = ph_data['col3'] #... define the time column

# photon times are turned into an array and converted to datetime format
np_times = np.array(ph_time)
timeincxo = Time(np_times, format='cxcsec')
chandra_evt_time = timeincxo.iso
# Chandra time then converted to a plotable format
chandra_evt_time = Time(chandra_evt_time, format='iso', out_subfmt='date_hm')
plot_time = Time.to_datetime(chandra_evt_time)
print('')
print('All observation will be analysed')

# Performing the coord transformation on the photons within the selected region

# The coordinate transformation is performed on the full observation

cxo_ints = []
sup_props_list = []
sup_time_props_list = []
sup_lat_list = []
sup_lon_list = []
lonj_max = []
latj_max = []
sup_psf_max = []
ph_tevts = []
ph_xevts = []
ph_yevts = []
ph_chavts = []
ph_sampvts = []; ph_sumampvts = []; ph_pivts = []; ph_ampsfvts = []
ph_av1vts = []; ph_av2vts = []; ph_av3vts = []
ph_au1vts = []; ph_au2vts = []; ph_au3vts = []
emiss_evts = []
ph_cmlevts = []
psfmax =[]


# perform the coordinate transformation for entire observation
tevents = ph_data['col3']
xevents = ph_data['col1']
yevents = ph_data['col2']
chaevents = ph_data['col4']
sampevents = ph_data['col5']; sumampsevents = ph_data['col6']; pievents = ph_data['col7']; ampsfevents = ph_data['col8']
av1events = ph_data['col9']; av2events = ph_data['col10']; av3events = ph_data['col11']
au1events = ph_data['col12']; au2events = ph_data['col13']; au3events = ph_data['col14']

"""CODING THE SIII COORD TRANSFORMATION"""
# define the local time and central meridian latitude (CML) during the observation  
jup_time = (eph_DOYFRAC_jup - evt_DOYFRAC)*86400.0 + tstart
jup_cml_0 = float(eph_jup['PDObsLon'][0]) + j_rotrate * (jup_time - jup_time[0])
interpfunc_cml = interpolate.interp1d(jup_time, jup_cml_0)
jup_cml = interpfunc_cml(tevents)
jup_cml = np.deg2rad(jup_cml % 360)
# find the distance between Jupiter and Chandra throughout the observation, convert to km
interpfunc_dist = interpolate.interp1d(jup_time, (eph_jup['delta'].astype(float))*AU_2_km)
jup_dist = interpfunc_dist(tevents)
dist = sum(jup_dist)/len(jup_dist)
kmtoarc = np.rad2deg(1.0/dist)*3.6E3 # convert from km to arc
kmtopixels = kmtoarc/scale # convert from km to pixels using defined scale
rad_eq_0 = 71492.0 # radius of equator in km
rad_pole_0 = 66854.0 # radius of poles in km
ecc = np.sqrt(1.0-(rad_pole_0/rad_eq_0)**2) # oblateness of Jupiter 
rad_eq = rad_eq_0 * kmtopixels
rad_pole = rad_pole_0 * kmtopixels # convert both radii form km -> pixels
alt0 = alt * kmtopixels # altitude at which we think emission occurs - agreed in Southampton Nov 15th 2017

# find sublat of Jupiter during each Chandra time interval
interpfunc_sublat = interpolate.interp1d(jup_time, (sub_obs_lat_jup.astype(float)))
jup_sublat = interpfunc_sublat(tevents)
# define the planetocentric S3 coordinates of Jupiter 
phi1 = np.deg2rad(sum(jup_sublat)/len(jup_sublat))
nn1 = rad_eq/np.sqrt(1.0 - (ecc*np.sin(phi1))**2)
p = dist/rad_eq
phig = phi1 - np.arcsin(nn1 * ecc**2 * np.sin(phi1)*np.cos(phi1)/p/rad_eq)
h = p * rad_eq *np.cos(phig)/np.cos(phi1) - nn1
interpfunc_nppa = interpolate.interp1d(jup_time, (eph_jup['NPole_ang'].astype(float)))
jup_nppa = interpfunc_nppa(tevents)
gamma = np.deg2rad(sum(jup_nppa)/len(jup_nppa))
omega = 0.0
Del = 1.0

#define latitude and longitude grid for entire surface
lat = np.zeros((int(360) // int(Del))*(int(180) // int(Del) + int(1)))
lng = np.zeros((int(360) // int(Del))*(int(180) // int(Del) + int(1)))
j = np.arange(int(180) // int(Del) + int(1)) * int(Del)

for i in range (int(0), int(360)):# // int(Del) - int(1)):
    lat[j * int(360) // int(Del) + i] = (j* int(Del) - int(90))
    lng[j * int(360) // int(Del) + i] = (i* int(Del) - int(0))

# perform coordinate transfromation from plentocentric -> planteographic (taking into account the oblateness of Jupiter
# when defining the surface features)
coord_transfo = gca_tools.ltln2xy(alt=alt0, re0=rad_eq_0, rp0=rad_pole_0, r=rad_eq, e=ecc, h=h, phi1=phi1, phig=phig, lambda0=0.0, p=p, d=dist, gamma=gamma,            omega=omega, latc=np.deg2rad(lat), lon=np.deg2rad(lng))

# Assign the corrected transformed position of the X-ray emission
xt = coord_transfo[0]
yt = coord_transfo[1]
cosc = coord_transfo[2]
condition = coord_transfo[3]
count = coord_transfo[4]

# Find latiutde and lonfitude of the surface features
laton = lat[condition] + 90
lngon = lng[condition]

# Define the limb of Jupiter, to ensure only auroral photons are selected for analysis
cosmu = gca_tools.findcosmu(rad_eq, rad_pole, phi1, np.deg2rad(lat), np.deg2rad(lng))
limb = np.where(abs(cosmu) < 0.05)

# This next step creates the parameters used to plot what is measured on Jupiter. In the code, I define this as "props" (properties)
# which has untis of counts/m^2. "timeprops" has units of seconds

# Creating 2D array of the properties and time properties
props = np.zeros((int(360) // int(Del), int(180) // int(Del) + int(1)))
timeprops = np.zeros((int(360) // int(Del), int(180) // int(Del) + int(1)))
num = len(tevents)
# define a Gaussian PSF for the instrument
psfn = np.pi*(fwhm / (2.0 * np.sqrt(np.log(2.0))))**2
# create a grid for the position of the properties
latx = np.zeros(num)
lonx = np.zeros(num)

# Equations for defining ellipse region
tilt_ang_rad = np.deg2rad(tilt_ang)
R_eq_as = (ang_diam/2.)/np.cos(tilt_ang_rad) # equatorial radius of Jupiter in arcsecs
R_pol_as = R_eq_as * np.sqrt(1 - ecc**2) # polar radius of Jupiter in arcsecs

for k in range(0,num-1):

    # convert (x,y) position to pixels
    xpi = (xevents[k]/scale)
    ypi = (yevents[k]/scale)

    # only select photons that lie within ellipse of Jupiter defined using JPL Horizons data
    if (xevents[k] * np.cos(tilt_ang_rad) + yevents[k] * np.sin(tilt_ang_rad)) ** 2./(R_eq_as ** 2) + (xevents[k] * np.sin(tilt_ang_rad) - yevents[k] * np.cos(tilt_ang_rad)) ** 2./(R_pol_as ** 2.) < 1.0:

        cmlpi = (np.rad2deg(jup_cml[k]))#.astype(int)

        xtj = xt[condition]
        ytj = yt[condition]
        latj = (laton.astype(int)) % 180
        lonj = ((lngon + cmlpi.astype(int) + 360.0).astype(int)) % 360
        dd = np.sqrt((xpi-xtj)**2 + (ypi-ytj)**2) * scale
        psfdd = np.exp(-(dd/ (fwhm / (2.0 * np.sqrt(np.log(2.0)))))**2) / psfn # define PSF of instrument

        psf_max_cond = np.where(psfdd == max(psfdd))[0] # finds the max PSF over each point in the grid
        count_mx = np.count_nonzero(psf_max_cond)
        if count_mx != 1: # ignore points where there are 2 cases of the same max PSF
            continue
        else:

            props[lonj,latj] = props[lonj,latj] + psfdd # assign the 2D PSF to the each point in the grid
            emiss = np.array(np.rad2deg(np.cos(cosc[condition[psf_max_cond]]))) # find the emission angle from each max PSF
            # record the corresponding photon data at each peak in the grid...
            emiss_evts.append(emiss)
            ph_cmlevts.append(cmlpi)
            ph_tevts.append(tevents[k])
            ph_xevts.append(xevents[k])
            ph_yevts.append(yevents[k])
            ph_chavts.append(chaevents[k])
            ph_sampvts.append(sampevents[k]); ph_sumampvts.append(sumampsevents[k]); ph_pivts.append(pievents[k]); ph_ampsfvts.append(ampsfevents[k])
            ph_av1vts.append(av1events[k]); ph_av2vts.append(av2events[k]); ph_av3vts.append(av3events[k])
            ph_au1vts.append(au1events[k]); ph_au2vts.append(au2events[k]); ph_au3vts.append(au3events[k])
            psfmax.append(psfdd[psf_max_cond])
            latj_max.append(latj[psf_max_cond])
            lonj_max.append(lonj[psf_max_cond])
            ph_tevts_arr = np.array(ph_tevts, dtype=float)
            ph_xevts_arr = np.array(ph_xevts, dtype=float)
            ph_yevts_arr = np.array(ph_yevts, dtype=float)
            ph_chavts_arr = np.array(ph_chavts, dtype=float)
            ph_sampvts_arr = np.array(ph_sampvts, dtype=float); ph_sumampvts_arr = np.array(ph_sumampvts, dtype=float); ph_pivts_arr = np.array(ph_pivts, dtype=float); ph_ampsfvts_arr = np.array(ph_ampsfvts, dtype=float)
            ph_av1vts_arr = np.array(ph_av1vts, dtype=float); ph_av2vts_arr = np.array(ph_av2vts, dtype=float); ph_av3vts_arr = np.array(ph_av3vts, dtype=float)
            ph_au1vts_arr = np.array(ph_au1vts, dtype=float); ph_au2vts_arr = np.array(ph_au2vts, dtype=float); ph_au3vts_arr = np.array(ph_au3vts, dtype=float)
            #... and save as text file
            np.savetxt(str(folder_path)+ "/%s_photonlist_full_obs_ellipse.txt" % obs_id, np.c_[ph_tevts_arr, ph_xevts_arr, ph_yevts_arr, ph_chavts_arr, latj_max, lonj_max, ph_cmlevts, emiss_evts, psfmax, ph_sampvts_arr, ph_sumampvts_arr, ph_pivts_arr, ph_ampsfvts_arr, ph_av1vts_arr, ph_av2vts_arr, ph_av3vts_arr, ph_au1vts_arr, ph_au2vts_arr, ph_au3vts_arr], delimiter=',', header="t(s),x(arcsec),y(arcsec),PHA,lat (deg),SIII_lon (deg),CML (deg),emiss (deg),Max PSF,samp,sumamps,pi,amp_sf,av1,av2,av3,au1,au2,au3", fmt='%s')

# effectively, do the same idea except for exposure time
obs_start_times = tevents.min()
obs_end_times = tevents.max()

interval = obs_end_times - obs_start_times

if interval > 1000.0:
    step = interval/100.0
elif interval > 100.0:
    step = interval/10.0
else:
    step = interval/2.0

time_vals = np.arange(round(int(interval/step)))*step + step/2 + obs_start_times

interpfunc_time_cml = interpolate.interp1d(jup_time,jup_cml_0)
time_cml = interpfunc_time_cml(time_vals)



for j in range(0, len(time_vals)):
    timeprops[((lngon + time_cml[j].astype(int))%360).astype(int),laton.astype(int)] = timeprops[((lngon + time_cml[j].astype(int))%360).astype(int),laton.astype(int)] + step


# record the fluxes and position of the max PSFs
sup_props_list = props
sup_time_props_list = timeprops

np.save(str(folder_path) + f'/{obsID}_sup_props_list.npy', np.array(sup_props_list))
np.save(str(folder_path) + f'/{obsID}_sup_time_props_list.npy', np.array(sup_time_props_list))


