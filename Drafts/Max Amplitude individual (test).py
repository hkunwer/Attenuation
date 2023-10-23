#!/usr/bin/env python
# coding: utf-8

# In[77]:


import os
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,arcsin,sqrt,abs,pi,log10,exp
from scipy.fftpack import fft,ifft

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
get_ipython().run_line_magic('matplotlib', 'inline')

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

from obspy.clients.neic import Client as nClient
from obspy.clients.fdsn import Client as fdsnClient

from obspy.geodetics.base import locations2degrees as l2d
from obspy.geodetics.base import degrees2kilometers as d2km
from obspy.geodetics.base import kilometers2degrees as km2d
from obspy.geodetics.base import gps2dist_azimuth as ll2az

from obspy.core.inventory.response import ResponseListElement as amplitude

#rtergpy
from rtergpy.waveforms import get_respinv, getwaves,process_waves, trstat2pd
from rtergpy.run import defaults, event

from locale import setlocale
from scipy.fftpack import fft,ifft
from scipy.stats import gmean 
from tqdm import tqdm
from compress_pickle import dump as cpkldump # reading/writing compressed pickles
from compress_pickle import load as cpklload # reading/writing compressed pickles


# In[3]:


Defaults=defaults()
Event=event()
Defaults.src="RASPISHAKE"
Defaults.network="AM"
Defaults.chan="EHZ"
Defaults.stationrange=[1.,10.]

Event.origin=[[18.4578,-73.3389,10],UTCDateTime(2022,1,24,13,16,23.425)]

st,df=getwaves(Defaults=Defaults,Event=Event)


# In[29]:


#invp = pd.read_csv('https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=2022-01-01T00:00:00&endtime=2022-12-29T00:00:00&latitude=17.5&longitude=-70&maxradius=5&minmag=4')
#print(invp)


# In[79]:


#Filtering using rtergpy process waves function

fmid=0.25
halfrange=0.25/2
fmin=fmid-halfrange
fmax=fmid+halfrange
stp=process_waves(st,freqmin=fmin,freqmax=fmax)
st[0].plot();
stp[0].plot();
print(stp[0])


# In[89]:


#Now try an pull max amplitudes


trdf=pd.DataFrame()
for tr in st:
        trdf=trdf.append(trstat2pd(tr),ignore_index=True)

trf=fft(trdf)
print(trdf)


# In[90]:


# from obspy 
# AmplitudeCategory = Enum([ "mean", "period", "integral"])
# AmplitudeUnit = Enum(["m", "s", "m/s", "m/(s*s)", "m*s", "dimensionless", "other"])


# In[91]:


#Using signal processing github

signal = tr

fft_spectrum = np.fft.rfft(signal)
fft_spectrum

#freq = np.fft.rfftfreq(signal.size, d=1/sampFreq) #cant use signal.size (doesnt work :c )

fft_spectrum_abs = np.abs(fft_spectrum)

plt.plot(freq, fft_spectrum_abs)
plt.xlabel("frequency, Hz")
plt.ylabel("Amplitude, units")
plt.show()


# In[ ]:


#Attenuation = Amplitude(Frequency) + geometric spread

