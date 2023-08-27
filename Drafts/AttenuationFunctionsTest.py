# %%
import csv
import numpy as np
from scipy.fftpack import fft,ifft

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
    fftdata = np.fft.rfft(data, n=samp, t=time)

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