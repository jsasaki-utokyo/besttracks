# besttracks

![tracks plot](https://raw.githubusercontent.com/miniufo/besttracks/master/pics/Global_TC_tracks.png)


## 1. Introduction
Tropical cyclone (TC) best-track datasets are analyzed, maintained, and hosted by several Regional Specialized Meteorological Centers (RSMCs), agencies, or projects all over the world.  These agencies include:
-  **JTWC:** Joint Typhoon Warning Center, Naval Oceanography Portal.  This agency currently hosted TC datasets over several ocean basins except the North Atlantic Ocean i.e.,  western Pacific basin (BWP), North Indian Ocean (BIO), and Southern Hemisphere basin (BSH).
https://www.metoc.navy.mil/jtwc/jtwc.html?best-tracks
-  **CMA:** China Meteorological Administration.  This agency only hosted the TC dataset over the western North Pacific.
http://tcdata.typhoon.org.cn/en/zjljsjj_zlhq.html
- **JMA:** RSMC Tokyo-Typhoon Center, Japan Meteorological Agency.  This agency only hosted the TC dataset over the western North Pacific
https://www.jma.go.jp/jma/jma-eng/jma-center/rsmc-hp-pub-eg/trackarchives.html
- **NHC:** National Hurricane Center, National Oceanic and Atmospheric Administration.  This agency hosted the TC datasets for both the North Atlantic Ocean and eastern North Pacific, which are not covered by JTWC.
https://www.nhc.noaa.gov/data/#hurdat
- **IBTrACS:** International Best Track Archive for Climate Stewardship.  This project merges the best-track datasets already exist at other agencies (more than the above) into a worldwide TC database.
https://www.ncdc.noaa.gov/ibtracs/


| RSMC | WNP | NEP | NAT | NIO | SIO | WSP | SAT |
| :----: | :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| JTWC | X |  |  | X | X | X |  |
| CMA | X |  |  |  |  |  |  |
| JMA | X |  |  |  |  |  |  |
| NHC |  | X | X |  |  |  |  |
| IBTrACS | X | X | X | X | X |  | X |


Unfortunately, different agencies use different data formats.  This python-based project **`besttracks`**, aims to provide a unified interface to access these datasets in different formats, and organizes them into a unified data structure called **`TCSet`** and **`TC`**, which are based on `pandas.DataFrame` that are convient for python users.  Simple plot of track and intensity is also easy and some basic statistics are also provided.

Up to now, the datasets from the above agencies are supported.  It would be good to add more agencies and more formats.  We also provide the parser function for CMA operational forecast data (BABJ format), which is also popular in China.

---

## 2. How to install on Windows 10
**Requirements**
`besttracks` is developed under the environment with `numpy` (=version 1.15.4), `pandas` (=version 1.0.3), `xarray` (=version 0.15.1), `matplotlib` (=version 3.3.1), and `cartopy` (=version 0.18.0).  Older versions of these packages are not well tested.

The following is based on [pyfvcom](https://github.com/jsasaki-utokyo/pyfvcom) environment; otherwise create a conda virtual environment, activate it, and install the prerequisite packages.

```
# Move into an appropriate folder where besttracks package is to be installed.
# The following is based on pyfvcom environment; otherwise activate a virtual environment and install prerequisite packages before git clone.
conda activate pyfvcom
conda install xarray
conda install scikit-learn
git clone https://github.com/jsasaki-utokyo/besttracks.git
cd besttracks
pip install -e .
```

---

## 3. Examples
### 3.0 Preparation
```python
from besttracks import parse_TCs
from besttracks import parseJMA

# Suppress deprecated warnings for cartopy 0.18.0
import warnings
warnings.simplefilter('ignore')

# Specify JMA best track data path
fpath = "D:/Github/besttracks/data/jma_rsmc/bst_all.txt"
```
### 3.1 Parsing JMA best-tracks into df (pandas.DataFrame)
`df.columns=['ID', 'NAME', 'TIME', 'LAT', 'LON', 'TYPE', 'PRS', 'WND']` where `ID`: 6 digits of year in 4-digit and typhoon number in 2-digit. `TIME`: in UTC, `WND` in knot (1 knot = 0.51444 m/s)

```python
df = parseJMA(fpath)
ID = '201919'
df[df.ID == ID].head()
```
### 3.2 Best-track datasets manipulations
Parsing best-track dataset **JMA** into `TCSet` would be as simple as:
```python
from besttracks import parse_TCs

# parse dataset from JMA
TCs_JMA = parse_TCs(fpath, agency='JMA')

# Brief describe the dataset
print(TCs_JMA)

# Plotting all TC tracks
TCs_JMA.plot_tracks(add_legend=True)
```

![tracks plot](https://raw.githubusercontent.com/miniufo/besttracks/master/pics/tracks_cma.png)

One can also bin the tracks into gridded statistics (also known as PDF distribution) as:
```python
# binning the tracks into gridded data
TCs_JMA.binning()
```

![binning plot](https://raw.githubusercontent.com/miniufo/besttracks/master/pics/binning_cma.png)

---

### 3.3 A single TC manipulation
Manipulating a single `TC` is also simple:
```python
# Selecting a single TC
tc = TCs_JMA[-1]

# Briefly descibe the TC
print(tc)

# Plot the TC track and intensity
tc.plot()
```
![tc plot](https://raw.githubusercontent.com/miniufo/besttracks/master/pics/tc_plot.png)

---

### 3.4 Timeseries statistics
`TCSet` also supports statistical analysis over time space. One can plot the timeseries of TC number and accumulated cyclonic energy (ACE) of a `TCSet` as:
```python
# plot the climatological timeseries of No. and ACE
TCs_JMA.plot_timeseries(freq='annual')
```
![tc plot](https://raw.githubusercontent.com/miniufo/besttracks/master/pics/timeseries.png)

More examples can be found at this [notebook](https://github.com/miniufo/besttracks/blob/master/notebooks/QuickGuide.ipynb)
