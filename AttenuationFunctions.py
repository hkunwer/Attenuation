import os
from locale import setlocale
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import abs
import pickle

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
    stp.remove_response(output="DISP") #using for comparison, will maybe use VEL for energies later. 
    stp.detrend(type='polynomial', order=5) # post-instrument removal
    stp.taper(taper)
    print("Completed: instrument response removal and taper, data returned as Displacement")
    
    return stp

def maxamp_calc(stream, EQ, Defaults, etime, eloc):
    
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
            threshold = 50 # Add a condition to exclude values above threshold
            if maxamp <= threshold:
                dist = tr.stats.distance / 1000
                station = tr.stats.station
                maxamps.append(maxamp)
                dist_str.append(dist)
                station_name.append(station)
                frequencies.append(freq_range)
  
    print("Completed: frequency band filtering")
          
    df_freq = pd.DataFrame({"maxamps":maxamps,"distance":dist_str,"station":station_name, "frequency":frequencies})         
            
    stations = []
    distances = []
    event_dates = []
    magnitudes = []
    magnitude_types = []
    max_amplitudes = []
    frequency_bands = []
    event_location = []

    for station, distance, max_amp, freq_range in zip(df_freq['station'], df_freq['distance'], df_freq['maxamps'], 
                                                      df_freq['frequency']):
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
        'magnitude': magnitudes,
        'mag type': magnitude_types,
        'max amplitude': max_amplitudes,
        'frequency band': frequency_bands
    }

    # Create a new DataFrame with the new data
    df_new = pd.DataFrame(data)
    print("Completed: made new dataframe of max amps for this event")
    
    return df_new

def update_and_save_dataframe(df_new):
    
    # Try to load an existing DataFrame if available, or create an empty one if not
    try:
        with open('AllResults.pkl', 'rb') as file:
            df_ALL = pickle.load(file)
            print("Loaded data from the AllResults.pkl file")
            
            # append df_new to the collective DataFrame
            df_ALL = df_ALL.append(df_new, ignore_index=False)

            # Save the updated DataFrame to a Pickle file
            try:
                with open('AllResults.pkl', 'wb') as file:
                    pickle.dump(df_ALL, file)
                    print("Finished updating AllResults.pkl:", df_ALL.tail())
            except Exception as e:
                print(f"Error while saving new info to Pickle file: {e}")
    
    except FileNotFoundError:
        df_ALL = pd.DataFrame(columns=["station name", "distance", "event", "location", "magnitude", "mag type", 
                                       "max amplitude", "frequency band"])
        print("Could not load data from existing file, creating new file named AllResults.pkl")

        # append df_new to the collective DataFrame
        
        df_ALL = df_ALL.append(df_new, ignore_index=False)
        # Save the updated DataFrame to a Pickle file
        try:
            with open('AllResults.pkl', 'wb') as file:
                pickle.dump(df_ALL, file)
                print("Finished updating AllResults_new.pkl:", df_ALL.tail())
        except Exception as e:
            print(f"Error while saving new info to Pickle file: {e}")
        
def maxamp_plot():
    
        with open('AllResults.pkl', 'rb') as pkl_file:
            dataframe = pickle.load(pkl_file)
            print("Completed: Loaded data from the AllResults.pkl file")

    color_mapping = {
    (0.10, 0.25): 'blue',
    (0.25, 0.50): 'green',
    (0.50, 0.75): 'red',
    (0.75, 1.00): 'orange',
    (1.00, 1.25): 'purple',}
    
    print("Plotting all max amps and distance now........")
          
    # Create a figure and axis
    fig, ax = plt.subplots()
    # Loop through unique frequency bands
    for freq_range, group in dataframe.groupby('frequency band'):
        color = color_mapping[tuple(freq_range)]  # Get color based on frequency range
        ax.scatter(group['distance'], group['max amplitude'], label = f'{freq_range[0]:.2f}-{freq_range[1]:.2f} Hz',
                   color=color, marker = 'o',facecolors='none')
    # Add labels, legend, and show the plot
    plt.yscale('log')
    plt.title('Max Amplitude for 2022-01-24')
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Max amplitude (Log scale)')
    plt.legend()
    image = plt.savefig('AllMaxAmps2022-01-24.png') 
    #save plot by date
    plt.show()  
    
    return dataframe