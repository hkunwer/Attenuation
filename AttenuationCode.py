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
from obspy import read
from obspy import read_inventory

from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")

from obspy.core.stream import Stream
from obspy.core.event import read_events
from obspy.core.inventory.inventory import read_inventory
from obspy.core.util import Enum
from obspy.core.inventory.response import ResponseListElement as amplitude
from obspy.core.event.magnitude import StationMagnitude

from obspy.clients.neic import Client as nClient
from obspy.clients.fdsn import Client as fdsnClient

from obspy.geodetics.base import locations2degrees as l2d
from obspy.geodetics.base import degrees2kilometers as d2km
from obspy.geodetics.base import kilometers2degrees as km2d
from obspy.geodetics.base import gps2dist_azimuth as ll2az

from obspy.signal.filter import bandpass

#rtergpy
from rtergpy.run import defaults, event, etime2name, src2ergs
from rtergpy.waveforms import get_respinv, getwaves, loadwaves, process_waves, trstat2pd

#attenuation
from AttenuationFunctions import processANSS, freqmaxes

# %%
# Processing and Reading information about event stored in ANSS_data.txt

Defaults=defaults()
Event=event()
Defaults.src='RASPISHAKE'
Defaults.network ="AM"
network = Defaults.network
Event.newData=False # True would collect new data. False if already downloaded.
Event.ecount='00'
edateold=""
processANSS() #process to remove unneccesary information
ANSS=pd.read_csv('ANSS_processed_data.csv', sep=',', comment="#")
for index, EQ in ANSS.iterrows():
    eloc = [EQ.LAT,EQ.LONG,EQ.DEPTH] 
    year,mo,dd = EQ.DATE.split('/')
    hh,mn,sec = EQ.TIME.split(':')
    etime=(UTCDateTime(int(year),int(mo),int(dd),int(hh),int(mn),float(sec)))
    # iterate ecount
    if EQ.DATE == edateold:
        Event.ecount=str(int(Event.ecount)+1).zfill(2)
    else:
        Event.ecount='00'
    edateold=EQ.DATE
    Event.eventname=etime2name(etime,ecount=Event.ecount)
    Event.origin=[eloc,etime]
    Event.focmech=[EQ.STK, EQ.DP, EQ.RKE] # phi,delta,lmbda

    print("\n\n"+Event.eventname+" ===============================")
    try:
        src2ergs(Defaults=Defaults,Event=Event)
    except:
        print("ERROR: running on "+Event.eventname+" failed!!!!\n\n")
print(ANSS) #to test if reading correctly

#somewhere here it stops working :( Am trying best to troubleshoot. 
#I think the issue is how the data is read, BUT it is printing ANSS correctly, so I think maybe the 
#actually assigning isnt working out? Will ask Andy

# %%
# getting waveforms is currently not running since switching to ANSS. 
Event.newData = True  # False means locally stored

if Event.newData:
    print("Getting waveforms")
    #st,df=getwaves(eloc,etime,pwindow=Def.waveparams[1],rads=Def.stationrange)
    st,df=getwaves(Defaults=Defaults,Event=Event)
else:
    print("Loading locally stored waveforms")
    try:
        st,df=loadwaves(Defaults=Defaults,Event=Event) 
    except:
        print("Couldn't load data for "+eventname+". Attempting to download:")
        st,df=getwaves(Defaults=Defaults,Event=Event)
if len(st) == 0:
    raise Exception("ERROR: No waveforms retreived.") 

# %%
# filter for instrument response and taper

taper=0.05

stp=st.copy()  # create backup
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
stp_freq5.filter("bandpass", freqmin=1, freqmax=1.25)

#add all freq to a tuple
stp_freqs = (stp_freq1,stp_freq2,stp_freq3,stp_freq4,stp_freq5)

# %%
# Plot all the wiggles, shows how frequency filtering works. 

stp[0].plot(); #only instrument response removal
stp_freq1[0].plot(); # bandpass filter of 0.1-0.25
stp_freq2[0].plot(); # bandpass filter of 0.25-0.5
stp_freq3[0].plot(); # bandpass filter of 0.5-0.75
stp_freq4[0].plot(); # bandpass filter of 0.75-0.1
stp_freq5[0].plot(); # bandpass filter of 1-1.25

# %%
# use freqmaxes function to get maxamps for all traces in one frequency band
# compare the plots formed from the function above with previous tests in Appendix below

ALLmaxamps, ALLidx_of_max_amps = [], []

for st in stp_freqs:

    for tr in st:
        
        # print(f"Now processing trace {tr}")
        maxamps, idx_of_max_amps = freqmaxes(tr, freqs)
        # print(maxamps, idx_of_max_amps) # to test if it worked
        ALLmaxamps.append(maxamps)
        ALLidx_of_max_amps.append(idx_of_max_amps)
        # trname.append(tr)  # no need to save this...it'll get huuuuuuuge
        # stname.append(st)
       
df = pd.DataFrame({"maxamps":ALLmaxamps, "idx of max amps":ALLidx_of_max_amps}) 
df

# Next we will plot all of the frequencies against each other. 
# Ask andy about help for next steps after this?