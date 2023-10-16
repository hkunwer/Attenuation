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
from Attenuation.AttenuationFunctions import processANSS, freqmaxes

# %%
# Processing and Reading information about event stored in ANSS_data.txt

Defaults = defaults()
Event = event()
Defaults.src="RASPISHAKE"
Defaults.network="AM"
Defaults.chan="EHZ"
Defaults.stationrange=[25.,30.]
Event.ecount='00'
Event.iter='RS'
#Event.newData=False   # use already downloaded data
Event.newData=True
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
        src2ergs(Defaults=Defaults,Event=Event) #should i use this or just getwaves. Ask andy. 
    except:
        print("ERROR: running on "+Event.eventname+" failed!!!!\n\n")