#!/usr/bin/env python
# coding: utf-8

# In[2]:


# VERY BASIC IDEA OF WHAT WE DOING. 
# You got this hiba! Code isnt THAT BAD!!!

import matplotlib.pyplot as plt
import numpy as np

# Sample data: Replace this with your own data
offsets = np.array([10, 20, 30, 40, 50])  # Example offset values (distance from each station to source)
max_amplitudes = np.array([0.5, 0.8, 1.2, 0.9, 1.5])  # Example max amplitudes
frequency_bands = ["0.1-0.25", "0.25-0.5", "0.5-0.75", "0.75-1", "1-1.25"]  # Example frequency bands

# Plotting
plt.figure(figsize=(10, 6))

# Plot individual frequency bands using different colors
for i, band in enumerate(frequency_bands):
    plt.scatter(offsets, max_amplitudes, marker='o', label=band, alpha=0.7)
    offsets += 2  # Adjust this value to spread out the bands on the x-axis

# Common labels and title
plt.xlabel("Offset / Distance")
plt.ylabel("Max Amplitude")
plt.title("Maximum Amplitude vs. Offset for Different Frequency Bands")

plt.legend()
plt.grid(True)
plt.show()


# In[3]:


# %%
import os
from locale import setlocale
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,arcsin,sqrt,abs,pi,log10,exp
from scipy.fftpack import fft,ifft
from scipy.io import wavfile
from scipy.stats import gmean 
from tqdm import tqdm
from compress_pickle import dump as cpkldump # reading/writing compressed pickles
from compress_pickle import load as cpklload # reading/writing compressed pickles

#obspy
from obspy import UTCDateTime

from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")

from obspy.geodetics.base import locations2degrees as l2d
from obspy.geodetics.base import degrees2kilometers as d2km
from obspy.geodetics.base import kilometers2degrees as km2d
from obspy.geodetics.base import gps2dist_azimuth as ll2az

#rtergpy
from rtergpy.run import defaults, event, etime2name, src2ergs
from rtergpy.waveforms import getwaves

#attenuation
from AttenuationFunctionsTest import processANSS, freqmaxes

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth
# %%
import csv
import numpy as np
from scipy.fftpack import fft,ifft


# In[ ]:


# Functions

def processANSS():

    # Read the data from the file
    with open('ANSS_data.txt', 'r') as file:
        lines = file.readlines()

    # Remove the header line
    header = lines[0].strip().split()[1:-3]  # Exclude the last two columns
    #comment = ['#']
    header = ['Date', 'Time'] + header
    lines = lines[1:]

    # Modify the data
    modified_lines = []
    for line in lines:
        columns = line.split()

        # Extracting the date and time from the first column
        first_column = columns[0]
        date, time = first_column.split('T')

        modified_line = [date, time] + columns[1:6]  # Exclude the last two columns
        modified_lines.append(modified_line)

    modified_lines.insert(0, header)

    # Write the modified data to a CSV file (Could use pickle instead for better information storage) .pkl
    #with open('ANSS_processed_data.pkl', 'wb') as file:
        #pkl.dump(modified_lines, file)

    with open('ANSS_processed_data.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(modified_lines)

# %%
def freqmaxes(trace, freqs): 
    
    data = trace.data                           # the actual data
    samp = data.shape[-1]                       # number of sample point
    delta = trace.stats.delta                  # sample spacing 
    time = trace.times()                       # time axis 

    # ---- compute fft for our data
    fftdata = np.fft.rfft(data, n=samp)

    # ---- Fdat has real and complex part, to compute power spectrum we take the
    fftdata_abs = np.abs(fftdata)

    # ---- compute x-axis in frequency-domain for plotting
    xf = np.linspace(0.0, 1.0/(2.0*delta), int(samp/2)+1)

    #make variables to store the values in
    maxamps = [] #max amplitudes here
    idx_of_max_amps = [] #position of max amps here

    for freq in freqs: 
        
        #argmin find specific position here
        idx_of_min_freq = (np.abs(xf-freq[0])).argmin()
        idx_of_max_freq = (np.abs(xf-freq[1])).argmin()
        
        fftdata_abs_slice = fftdata_abs[idx_of_min_freq:idx_of_max_freq]

        # ---- get the max value of amplitude spectrum
        maxamp = fftdata_abs_slice.max()
        idx_of_max_amp = int(np.where(fftdata_abs_slice == maxamp) + idx_of_min_freq)
        maxamps.append(1.0/samp * maxamp)
        idx_of_max_amps.append(idx_of_max_amp)
    
    return maxamps, idx_of_max_amps

def get_respinv(network,eloc,etime,rads,chan,src="IRIS",**kwargs):
    """
    Algorithm to use Obspy's metadata to pull response and other metadata.  Returns an apporpiate
    inventory class.
    A.V. Newman Mon Jul 26 15:26:35 EDT 2021
    """
    client = fdsnClient(src)

    elat,elon,edep = eloc
    minrad,maxrad = rads 
    
    inventory = client.get_stations(network = network, 
            latitude = elat, longitude = elon, 
            minradius = minrad, maxradius = maxrad, 
            starttime=etime-86400, endtime=etime,
            channel=chan,
            #location="00",
            #matchtimeseries=True, # may not work if client is not iris
            #filename="test.txt", format="text",  # cannot be used with response
            #filename="test.xml", format="xml",
            level="response"
            )
    return inventory  # this is an Obspy type


# In[ ]:


Defaults = defaults()
Event = event()
Defaults.src="RASPISHAKE"
Defaults.network="AM"
Defaults.chan="EHZ"
Defaults.stationrange=[1.,10.]
Event.ecount='00'
Event.iter='RS'
# Event.newData = False   # use already downloaded data
Event.newData=False
edateold=""


# In[ ]:


# %%
# Processing and Reading information about event stored in ANSS_data.txt

#processANSS() #process to remove unneccesary information
# ANSS = pd.read_csv('ANSS_processed_data.csv', sep=',', comment='#')
ANSS = pd.read_csv('ANSS_processed_data.csv')
#print(ANSS) #Just to check if processing correctly

# run everything above to test on command line
for index, EQ in ANSS.iterrows():
    network = "AM"
    chan = "EHZ"
    src = "RASPISHAKE"
    rads = [1.,10.]
    eloc = [EQ.Latitude,EQ.Longitude,EQ.Depth] 
    MagType = [EQ.Mtype]
    MagValue = [EQ.Mag]
    Magnitude = [MagType, MagValue]
    year,mo,dy = EQ.Date.split('-')
    hh,mn,sec = EQ.Time.split(':')
    etime=(UTCDateTime(int(year),int(mo),int(dy),int(hh),int(mn),float(sec)))
    
    if EQ.Date == edateold:
        Event.ecount=str(int(Event.ecount)+1).zfill(2)
    else:
        Event.ecount='00'
    edateold=EQ.Date
    Event.eventname=etime2name(etime,ecount=Event.ecount)
    Event.origin=[eloc,etime]

    print("\n\n"+Event.eventname+" ===============================")
    try:
        st, df = [], []
        st, df = getwaves(Defaults=Defaults,Event=Event)
    except:
        print("ERROR: running on "+Event.eventname+" failed!!!!\n\n")
               
    inventory = get_respinv(network,eloc,etime,rads,chan,src)
    print(inventory)
    
    # Extract the latitude and longitude from the inventory
    station_latitude = inventory[0][0].latitude
    station_longitude = inventory[0][0].longitude
    
    # Now you have the latitude and longitude of the station
    print("Station Latitude:", station_latitude)
    print("Station Longitude:", station_longitude)

for index, EQ in ANSS.iterrows():
    # Event and station coordinates
    event_coords = (EQ.Latitude, EQ.Longitude)
    station_coords = (station_latitude, station_longitude)  # Replace with actual station coordinates

    # Calculate distance and azimuth
    distance, _, _ = gps2dist_azimuth(*event_coords, *station_coords)

    # Store the calculated distance in the ANSS DataFrame
    ANSS.at[index, "Distance"] = distance       


# In[ ]:


# filter for instrument response and taper

taper=0.05

stp = st.copy()  # create backup
# process data
stp.detrend(type='polynomial', order=5) # pre-instrument removal
stp.taper(taper)
stp.remove_response(output="DISP")
stp.detrend(type='polynomial', order=5) # post-instrument removal
stp.taper(taper)

# %%
# process each individual frequency band in a bandpass using .filter
# frequency bands of 0.1-0.25, 0.25-0.5, 0.5-0.75, 0.75-1, 1-1.5

freqs = ([0.1,0.25],[0.25,0.5],[0.50,0.75],[0.75,1],[1,1.25])

# freq1 band = 0.1-0.25,
stp_freq1 = stp.copy()
stp_freq1.filter("bandpass",freqmin=0.10 , freqmax=0.25)

# freq2 band = 0.25-0.5,
stp_freq2 = stp.copy()
stp_freq2.filter("bandpass", freqmin=0.25, freqmax=0.50)

# freq3 band = 0.5-0.75,
stp_freq3 = stp.copy()
stp_freq3.filter("bandpass", freqmin=0.50, freqmax=0.75)

# freq4 band = 0.75-1.
stp_freq4 = stp.copy()
stp_freq4.filter("bandpass", freqmin=0.75, freqmax=1.0)

# freq5 band = 1-1.25
stp_freq5 = stp.copy()
stp_freq5.filter("bandpass", freqmin=1.0, freqmax=1.25)

#add all freq to a tuple
stp_freqs = (stp_freq1,stp_freq2,stp_freq3,stp_freq4,stp_freq5)


# In[ ]:


# get max amps and offset

# use freqmaxes function to get maxamps for all traces in one frequency band
all_max_amplitudes = []

# Loop through each stream
for st in stp_freqs:

    maxamps_list = []  # Initialize empty list for max amplitudes of this stream
    idx_of_max_amps_list = []  # Initialize empty list for indices of max amplitudes of this stream

    # Loop through each trace within the stream
    for tr in st:
        maxamps, idx_of_max_amps = freqmaxes(tr, freqs)  # Calculate max amplitudes, discard indices
        maxamps_list.extend(maxamps)  # Extend max amplitudes list for this stream
    
    all_max_amplitudes.extend(maxamps_list)  # Extend max amplitudes list for all streams

# Create a DataFrame from the combined lists
df_amps = pd.DataFrame({"Max Amplitude": all_max_amplitudes})


# In[ ]:


# plot 

# Calculate offsets (assuming they are stored in ANSS DataFrame)
offsets = ANSS["Distance"]  # Adjust column name if needed

# Plotting
plt.figure(figsize=(10, 6))

# Loop through each frequency band and plot max amplitudes against offsets
for freq_band_index, (st, freq_range) in enumerate(zip(stp_freqs, freqs)):
    max_amplitudes_list = []  # Initialize empty list for max amplitudes
    
    # Loop through each trace within the frequency band
    for tr in st:
        maxamps, _, _ = freqmaxes(tr, [freq_range])  # Calculate max amplitudes for the specific frequency band
        max_amplitudes_list.extend(maxamps)  # Extend max amplitudes list for this trace
    
    # Scatter plot for max amplitudes against offsets
    plt.scatter(offsets, max_amplitudes_list, marker='o', label=f"Freq Band {freq_band_index + 1}")

# Common labels and title
plt.xlabel("Offset / Distance")
plt.ylabel("Max Amplitude")
plt.title("Maximum Amplitude vs. Offset for Different Frequency Bands")

plt.legend()
plt.grid(True)
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




