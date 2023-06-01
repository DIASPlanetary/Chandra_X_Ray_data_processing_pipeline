# Chandra_X_Ray_data_processing_pipeline
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5657141.svg)](https://doi.org/10.5281/zenodo.5657141)

This package contains the python files needed in the post-processing of data taken with the High Resolution Camera (HRC-I) instrument on-board the Chandra X-ray Observatory (CXO). 

The code is set up to process observations of Jupiter, but can be applied to other solar system objects in an analogous way. A description of each step in the pipeline is provided below, along with information on how to download the data.

## Contents
- [1. Downloading Data](#downloading-data)
   * [1.1. Pre-Processing](#pre-processing)
   * [1.2 chandra_epro](#chandra_repro)
- [2. SSO_FREEZE](#sso_freeze)
- [3. GO_CHANDRA](#go_chandra)      
- [4. PI_FILTER](#pi_filter)
- [5. Requirements](#requirements)


## Downloading Data

Chandra HRC-I data can be obtained directly from the [Chandra Data Archive](https://cda.harvard.edu/chaser/). 

Alternatively, there is a software package developed by the Chandra X-ray Center called CIAO, the Chandra Interactive Analysis of Observations, which has the capability of downloading the data directly from the command line. Instructions to download CIAO can be found here: https://cxc.cfa.harvard.edu/ciao/download/. 

> Note: the easiest way to download CIAO is through the Python distribution [Anaconda](https://www.anaconda.com/download).

Once CIAO has been installed, data from a specific observation can be downloaded using the ```download_chandra_obsid``` command (see https://cxc.cfa.harvard.edu/ciao/ahelp/download_chandra_obsid.html for documentation).

For example, data from the observation of Jupiter with the **observation ID** (obsID) 1862 can be dowloaded with the following command:

```shell
download_chandra_obsid 1862
```

### Pre-processing

Once downloaded, the data will be stored as **1862/**, with the sub-directories **primary/** and **secondary/**. The files needed for post-processing must be unpacked first using the following command:

```shell
gunzip primary/*.gz secondary/*.gz
```

### chandra_repro

The processing pipeline requires **header information** from the raw **.fits** files that is not present in some of the earlier observations of Jupiter. The obsIDs for these observations are listed in **no_samp.txt**. This step can be ignored for any obsIDs not present in this list.

Another CIAO command, [`chandra_repro`](https://cxc.cfa.harvard.edu/ciao/ahelp/chandra_repro.html), must be executed on the observations in **no_samp.txt** to reprocess the event files. When running `chandra_repro` the user will be propmted to provide input and output directories. These should be your data directory (e.g. .../1862/) and a new sub-directory (e.g. ...1862/repro/) respectively.

After performing ```chandra_repro```, this new third sub-directory, **repro/**, will contain the reprocessed event file. The orbital ephemeris file for the observation also needs to be moved into the **repro/** directory before running SSO_FREEZE. 

```shell
mv primary/orbit*_eph1.fits repro/
``` 

## SSO_FREEZE
The raw data obtained from HRC-I first have to be transformed into a frame of reference centered on Jupiter. The SSO_FREEZE algorithm uses appropriate ephemerides data from the JPL Horizons program and Chandra orbit ancillary data from the Chandra X-ray Center to account for Jupiter's motion on the sky and the relative position of the detector. The raw data are reprojected from sky x and y co-ordinates to a reference frame which is fixed to the motion of Jupiter.

```shell
python python_sso_freeze_v6.py
``` 

**Input:** uncorrected event file (hrcf*_evt2.fits)

**Output:** Jupiter centered event file (hrcf*_pytest_evt2.fits)

Note: Need to change the observation ID (**obsID**) and the **folder_path** variables in the **config.ini** file.

<sub>**Authors:** Dale Weigt, adapted from Randy Gladstone's 'stop_hrci' IDL script. Other contributors to the current (and/or) previous iteration(s) of the code are: Hunter Waite, Kurt Franke, Peter Ford, Seán McEntee, Caitríona Jackman, Will Dunn, Brad Snios, Ron Elsner, Ralph Kraft, and Graziella Branduardi-Raymont.</sub>

## GO_CHANDRA
Python script that takes the corrected .fits file from SSO_FREEZE and performs a coordinate transformation on the X-ray emission to wrap the point-spread function (PSF) around Jupiter. JPL Horizons data are used to define an ellipse which constrains the limb of the planet, enabling the selection of photons based on whether they lie on Jupiter's surface. The output of the GO_CHANDRA algorithm is a text file containing a list of all the time-tagged photons that emanate from Jupiter's surface. Information about the location of each photon (lat, SIII lon, CML) is contained within this text file, along with the amplifier signals needed to calculate the pulse invariant (PI) channels.

```shell
python python_go_chandra_v5.py
```

**Input:** Jupiter centered event file (hrcf*_pytest_evt2.fits)

**Output:** Time-tagged list of all jovian photons (*_photonlist_full_obs_ellipse.txt)

<sub>**Authors:** Dale Weigt, adapted from Randy Gladstone's 'go_chandra' IDL script. Other contributors to the current (and/or) previous iteration(s) of the code are: Hunter Waite, Kurt Franke, Peter Ford, Seán McEntee, Caitríona Jackman, Will Dunn, Brad Snios, Ron Elsner, Ralph Kraft, and Graziella Branduardi-Raymont.</sub>

## PI_FILTER
The gain on HRC-I has been degrading over time. To combat this, a new gain metric, pulse invariant (PI), was introduced which shifts the distribution of the amplifier signals to what was observed at the beginning of the Chandra mission. This step is important if the user wishes to compare observations over a long time-span. An attempt is made to minimise the particle background striking the detector by only including PI channels where the source (Jupiter) is expected to dominate the background. For this reason, only photons with PI values in the range 10-250 are included after this PI filter is applied (see https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022JA030971 for more detail).

```shell
python PI_filter_v2.py
```

**Input:** Time-tagged list of all jovian photons (*_photonlist_full_obs_ellipse.txt)

**Output:** Time-tagged list of jovian photons in PI range 10-250 (*_photonlist_PI_filter_Jup_full_10_250.txt)

<sub>**Authors:** Seán McEntee, Vinay Kashyap, Dale Weigt, Caitriona Jackman</sub>

## Requirements
* numpy 1.21.6
* pandas 1.2.3
* os 
* astropy 5.1
* astroquery 0.4.1
* scipy 1.7.3
* configparser
* matplotlib 3.5.2
* datetime
