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
from Drafts.AttenuationFunctionsTest import processANSS, freqmaxes

# %%
# Processing and Reading information about event stored in ANSS_data.txt

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

processANSS() #process to remove unneccesary information
# ANSS = pd.read_csv('ANSS_processed_data.csv', sep=',', comment='#')
ANSS = pd.read_csv('ANSS_processed_data.csv')
#print(ANSS) #Just to check if processing correctly

# run everything above to test on command line
for index, EQ in ANSS.iterrows():
    eloc = [EQ.Latitude,EQ.Longitude,EQ.Depth] 
    MagType = [EQ.Mtype]
    MagValue = [EQ.Mag]
    Magnitude = [MagType, MagValue]
    #These are just to test if reading correctly! 
    #print(MagType)
    #print(MagValue)
    #print("Magnitude:", Magnitude)
    #print(eloc)
    year,mo,dy = EQ.Date.split('-')
    hh,mn,sec = EQ.Time.split(':')
    etime=(UTCDateTime(int(year),int(mo),int(dy),int(hh),int(mn),float(sec)))
    # print("Time of event:", etime)
    # iterate ecount
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
# %%
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

# %%
# Plot all the wiggles, shows how frequency filtering works. 

#stp[0].plot(); #only instrument response removal
#stp_freq1[0].plot(); # bandpass filter of 0.1-0.25
#stp_freq2[0].plot(); # bandpass filter of 0.25-0.5
#stp_freq3[0].plot(); # bandpass filter of 0.5-0.75
#stp_freq4[0].plot(); # bandpass filter of 0.75-0.1
#stp_freq5[0].plot(); # bandpass filter of 1-1.25

# %%
# use freqmaxes function to get maxamps for all traces in one frequency band

# Initialize lists to store center frequencies, maximum amplitudes, and indices for each stream
all_center_frequencies = []
all_max_amplitudes = []

# Loop through each stream
for st in stp_freqs:

    maxamps_list = []  # Initialize empty list for max amplitudes of this stream
    idx_of_max_amps_list = []  # Initialize empty list for indices of max amplitudes of this stream
    center_frequencies_list = []  # Initialize empty list for center frequencies of this stream

    # Loop through each trace within the stream
    for tr in st:
        maxamps, idx_of_max_amps, center_frequencies = freqmaxes(tr, freqs)  # Calculate max amplitudes, discard indices
        maxamps_list.extend(maxamps)  # Extend max amplitudes list for this stream
        center_frequencies_list.extend(center_frequencies)  # Extend center frequencies list for this stream
    
    all_max_amplitudes.extend(maxamps_list)  # Extend max amplitudes list for all streams

# Create a DataFrame from the combined lists
df_amps = pd.DataFrame({"Center Frequency": all_center_frequencies, "Max Amplitude": all_max_amplitudes})

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


