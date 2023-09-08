import os
from locale import setlocale
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import abs

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

    ANSS = pd.read_csv('ANSS_processed_data.csv')
    
    return ANSS
        
def filtering(stream):
    
    taper=0.05
    stp = stream.copy()  # create backup
    # process data (maybe make this into a function)
    stp.detrend(type='polynomial', order=5) # pre-instrument removal
    stp.taper(taper)
    stp.remove_response(output="DISP")
    stp.detrend(type='polynomial', order=5) # post-instrument removal
    stp.taper(taper)
    
    return stp

def maxamp_calc(stream, eventID):
    
    freqlist = [(0.10, 0.25),(0.25, 0.50),(0.50, 0.75),(0.75, 1.00),(1.00, 1.25)]
    maxamps = []
    dist_str = []
    station_name = []
    frequencies = []

    for freq_range in freqlist:
        
        frequencymin, frequencymax = freq_range
        
        stp_freq = stream.copy()
        stp_freq.filter("bandpass", freqmin=frequencymin, freqmax=frequencymax)

        for tr in stp_freq:
            maxamp = np.max(abs(tr))
            dist = tr.stats.distance / 1000
            station = tr.stats.station
            maxamps.append(maxamp)
            dist_str.append(dist)
            station_name.append(station)
            frequencies.append(freq_range)
  
    df_freq = pd.DataFrame({"maxamps":maxamps,"distance":dist_str,"station":station_name, "frequency":frequencies}) 
    df_freq.to_csv('MaxAmps_'+eventID+'.csv', index=False)
    return df_freq

def organize_data(df_freq, EQ, etime,eloc, eventID):          
            
    stations = []
    distances = []
    event_dates = []
    magnitudes = []
    magnitude_types = []
    max_amplitudes = []
    frequency_bands = []
    event_location = []

    for station, distance, max_amp, freq_range in zip(df_freq['station'], df_freq['distance'], df_freq['maxamps'], df_freq['frequency']):
        stations.append(station)
        distances.append(distance)
        event_dates.append(etime)
        event_location.append(eloc)
        magnitudes.append(EQ.Mag)
        magnitude_types.append(EQ.Mtype)
        max_amplitudes.append(max_amp)
        frequency_bands.append(freq_range)

    data = {
        'station name': stations,
        'distance': distances,
        'event': event_dates,
        'location':event_location,
        'magnitude of event': magnitudes,
        'type of magnitude': magnitude_types,
        'max amplitude': max_amplitudes,
        'frequency band': frequency_bands
    }
    df_final = pd.DataFrame(data)
    df_final.to_csv('Results_'+eventID+'.csv', index=False) # Export to CSV

def maxamp_plot(eventID, dataframe):
    color_mapping = {
    (0.10, 0.25): 'blue',
    (0.25, 0.50): 'green',
    (0.50, 0.75): 'red',
    (0.75, 1.00): 'orange',
    (1.00, 1.25): 'purple',}
    # Create a figure and axis
    fig, ax = plt.subplots()
    # Loop through unique frequency bands
    for freq_range, group in dataframe.groupby('frequency'):
        color = color_mapping[tuple(freq_range)]  # Get color based on frequency range
        ax.scatter(group['distance'], group['maxamps'], label = f'{freq_range[0]:.2f}-{freq_range[1]:.2f} Hz',color=color, marker = 'o',facecolors='none')
    # Add labels, legend, and show the plot
    plt.yscale('log')
    plt.title('The Max Amplitude for '+ eventID)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Max amplitude (Log scale)')
    plt.legend()
    image = plt.savefig('MaxAmps_'+eventID+'.png') #save plot by unique event ID
    #plt.show()  