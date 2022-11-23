# File that takes the output photonlist of go_chandra as an input and filters based on the Pulse Invariant (PI). Only photons having correspnoding PI values in the range 10-250 are selected. This aims to account for the gain degradation over time in the High Resolution Camera (HRC) on board the Chandra X-ray Observatory (CXO).

# Authors: Seán McEntee, Vinay Kashyap, Dale Weigt, Caitríona Jackman.

#relevant packages 
import numpy as np
import pandas as pd
from astropy.time import Time

# Reading in configuration file containing any hard wired inputs.
import configparser
config = configparser.ConfigParser()
config.read('config.ini')

# Select observation ID (obsID)
obsID = config['inputs']['obsID']

# Accounting for different filepaths of ObsIDs that originlly had SAMP values and others that did not.
df = pd.read_csv('ObsIDs_with_samp.txt', header=None, delimiter='\t')
samp_ids = np.array(df.iloc[:,0])

# if int(obsID) in samp_ids:
#     folder_path = '/Users/mcentees/Desktop/Chandra/' + str(obsID) + '/primary'
# else:
#     folder_path = '/Users/mcentees/Desktop/Chandra/' + str(obsID) + '/repro'

folder_path = config['inputs']['folder_path']

# Reading in ellipse data
ellipse_data = pd.read_csv(folder_path + f'/{obsID}_photonlist_full_obs_ellipse.txt')

# reading in amplifier signals
av1_jup = np.array(ellipse_data['av1'])
av2_jup = np.array(ellipse_data['av2'])
av3_jup = np.array(ellipse_data['av3'])

au1_jup = np.array(ellipse_data['au1'])
au2_jup = np.array(ellipse_data['au2'])
au3_jup = np.array(ellipse_data['au3'])

# reading in amplifier scale factor 
amp_sf_jup = np.array(ellipse_data['amp_sf'])

# calculating sumamp values
sumamp_jup = av1_jup + av2_jup + av3_jup + au1_jup + au2_jup + au3_jup

# calculating samp values
samp_jup = (sumamp_jup * (2. ** (amp_sf_jup - 1.0)))/148.0

# Decimal year of beginning of observation needed to perform PI calculation - can be obtained from reading in catalogue containing key info for all observations.
# catalogue_path = config['PI Filter']['catalogue_path']
chandra_props = pd.read_excel('catalogue_all_data.xlsx')
index = np.where(chandra_props['ObsID'] == int(obsID))[0][0]
date_start = chandra_props['Start Date'][index]
date_dec = Time(date_start).decimalyear

# PI calculation
g= 1.0418475 + 0.020125799 * (date_dec - 2000.) + 0.010877227 * (date_dec - 2000.) ** 2. + - 0.0014310146 * (date_dec - 2000.) ** 3. + 5.8426766e-05 * (date_dec - 2000.) ** 4. # multiplicative scaling factor - Provided by Vinay Kashyap from the Chandra calibration team.
PI_jup = g * samp_jup

# creating output .txt file
df_filtered = ellipse_data.drop(['PHA','samp', 'sumamps', 'pi'], axis=1) # removing quantities that have been recalculated
df_filtered.insert(15, "PI", PI_jup) # adding the recalculated PI values to a new column
df_filtered['lat (deg)'] = df_filtered['lat (deg)'] - 90. # change latitude range to [-90,90]
df_filtered_Jup_full = df_filtered[(df_filtered["PI"] > 10) & (df_filtered["PI"] < 250)] # only photons that lie in the PI channel range 10-250 are selected.
df_filtered_Jup_full.to_csv(folder_path + f'/{obsID}_photonlist_PI_filter_Jup_full_10_250.txt', index=False) # writing PI filtered  photon list to new file

