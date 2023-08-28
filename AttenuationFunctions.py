# %%
import csv
import numpy as np
from scipy.fftpack import fft,ifft
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

#rtergpy
from rtergpy.run import defaults, event, etime2name, src2ergs
from rtergpy.waveforms import getwaves, get_respinv

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
    time = trace.times()       
                # time axis 
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
# need to fix this function so it actually works when called. Currently doing manually since 
# shes rude right now
def maximumAmplitude(stream, frequencymin, frequencymax):
    
    stp_freq = stream.copy()
    stp_freq.filter("bandpass",freqmin=frequencymin , freqmax= frequencymax)

    # stp_freq[0].plot();

    maxamps = [] 
    dist_str = []

    for tr in stp_freq:

        maxamp = np.max(abs(tr))
        dist = tr.stats.distance / 1000
        maxamps.append(maxamp)
        dist_str.append(dist)
        
    return maxamps, dist_str