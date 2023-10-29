import os
from locale import setlocale
import numpy as np
import csv
from numpy import abs
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

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
    
    try:
        # Try to load an existing DataFrame if available
        with open('Results.pkl', 'rb') as file:
            df_ALL = pickle.load(file)
            print("Loaded data from the Results.pkl file")
    except (FileNotFoundError, EOFError):
        # If the file doesn't exist or is not valid, create an empty DataFrame
        df_ALL = pd.DataFrame(columns=["station name", "distance", "event", "location", "magnitude", "mag type", "max amplitude", "frequency band"])
        print("Could not load data from an existing file, creating a new file named Results.pkl")

    # Append df_new to the collective DataFrame
    df_ALL = df_ALL.append(df_new, ignore_index=True)

    # Save the updated DataFrame to a Pickle file
    with open('Results.pkl', 'wb') as file:
        pickle.dump(df_ALL, file)
        print("Finished updating Results.pkl")
        
def normalization(): #removes outliers and normalizes data
    
    try: 
    # Load the original DataFrame from the pickle file
        with open('Results.pkl', 'rb') as pkl_file:
            original_dataframe = pickle.load(pkl_file)
            print("Successfully read Results.pkl")

    # outlier removal
        print("Removing outliers")
        # Calculate Z-scores for log-transformed max amplitude data
        z_scores = np.abs(stats.zscore(np.log(original_dataframe['max amplitude'])))

        # Set a threshold for outlier detection (e.g., Z-score greater than 3)
        threshold = 3

        # Remove outliers based on the threshold
        outlier_removed_df = original_dataframe[(z_scores < threshold)]

        # Filter only the frequency band (1.0, 1.25) as a reference
        filtered_dataframe = outlier_removed_df[outlier_removed_df['frequency band'] == (1.0, 1.25)]

    # normalization factor
        print("Beginning normalization")
        # Calculate the best fit line for the log-transformed filtered data
        filtered_distance = filtered_dataframe['distance']
        filtered_max_amplitude = filtered_dataframe['max amplitude']
        filtered_log_max_amplitude = np.log10(filtered_dataframe['max amplitude'])
        filtered_coefficients = np.polyfit(filtered_distance, filtered_log_max_amplitude, 1)
        filtered_best_fit_line = np.poly1d(filtered_coefficients)

        # Calculate the normalization factor (max amps) at a distance of 500 using the best-fit line
        normalizeBench = 100
        normalization_factor = filtered_best_fit_line(normalizeBench)
        #print(filtered_best_fit_line(500))
        
        print("creating normalized dataframe")

        # Create a new DataFrame for the normalized data
        normalized_dataframe = outlier_removed_df.copy()

        def normalize_max_amplitude(row, normalization_factor):
            return np.log10(row['max amplitude']) - normalization_factor 
        
        # Apply the normalization function to create the 'max amplitude (normalized)' column
        normalized_dataframe['max amplitude (normalized)'] = normalized_dataframe.apply(lambda row: normalize_max_amplitude(row, normalization_factor), axis=1)

        print("saved normalized.pkl")
        # Save the normalized DataFrame to a pkl file
        normalized_dataframe.to_pickle('normalized.pkl')
    #Plotting
        print("Plotting figures now")
        # Plot the log-transformed original data
        fig, axs = plt.subplots(1, 3, figsize=(18, 4))  # Use sharey=True
        #axs[1].sharey(axs[2])

        # Subplot 1: Filtered data with best fit line
        star_x = normalizeBench
        star_y = filtered_best_fit_line(star_x)
        axs[0].scatter(star_x, star_y, color='red', marker='*', s=50, label='Normalization Factor')

        axs[0].scatter(filtered_distance, filtered_log_max_amplitude, marker='o', label='(1.0,1.25)Hz', alpha=0.5)
        axs[0].plot(filtered_distance, filtered_best_fit_line(filtered_distance), color='red', label='Best Fit Line')
        axs[0].set_title('Filtered Data')
        axs[0].set_xlabel('Distance')
        axs[0].set_ylabel('Log10(Max Amplitude)')

        # Add equation of the line as a text annotation (adjust position)
        equation = f'Line Equation: y = {filtered_best_fit_line[0]:.3f}x + {filtered_best_fit_line[1]:.3f}'
        axs[0].text(0.30, 0.10, equation, transform=axs[0].transAxes, fontsize=9)  # Adjust position here
        axs[0].legend()

        # Subplot 2: Original data with different colors for each frequency band (excluding 1.0,1.25 Hz)
        color_mapping = {
            (0.10, 0.25): 'blue',
            (0.25, 0.50): 'green',
            (0.50, 0.75): 'red',
            (0.75, 1.00): 'orange',
            (1.00, 1.25): 'purple',
        }
        for freq_range, group in normalized_dataframe.groupby('frequency band'):
            color = color_mapping[tuple(freq_range)]
            axs[1].scatter(group['distance'], np.log10(group['max amplitude']), marker='o', facecolors='none',
                           alpha=0.5, label=f'{freq_range[0]:.2f}-{freq_range[1]:.2f} Hz', color=color)
        axs[1].set_title('Original Data')
        axs[1].set_xlabel('Distance')
        axs[1].set_ylabel('Log10(Max Amplitude)')
        axs[1].legend()

        # Subplot 3: Normalized data with different colors for each frequency band (including 1.0,1.25 Hz)

        star_x = normalizeBench
        star_y = filtered_best_fit_line(star_x) - filtered_best_fit_line(star_x)
        axs[2].scatter(star_x, star_y, color='red', marker='*', s=50, label='Normalization Factor')

        for freq_range, group in normalized_dataframe.groupby('frequency band'):
            color = color_mapping[tuple(freq_range)]
            axs[2].scatter(group['distance'], group['max amplitude (normalized)'], marker='o', facecolors='none',
                           alpha=0.5, label=f'{freq_range[0]:.2f}-{freq_range[1]:.2f} Hz', color=color)
        axs[2].set_title('Normalized Data')
        axs[2].set_xlabel('Distance')
        axs[2].set_ylabel('Log10(Max Amplitude) (Normalized)')
        axs[2].legend()
        print("saving FinalPlot.png")
        #image3 = plt.savefig("FinalPlot.png")
        image_filename = os.path.join('plots', 'FinalPlot.png')
        plt.savefig(image_filename)
        plt.tight_layout()
        plt.show()

        # Create the "plots" folder if it doesn't exist
        if not os.path.exists('plots'):
            os.mkdir('plots')

        return normalized_dataframe
    
    except: 
        print("Could not read Results.pkl or plot figure")
        return original_dataframe

def maxamp_plot():
    
    # Create the "plots" folder if it doesn't exist
    if not os.path.exists('plots'):
        os.mkdir('plots')
    
    try: 
        with open('Results.pkl', 'rb') as pkl_file:
            dataframe = pickle.load(pkl_file)
            print("Completed: Loaded data from the Results.pkl file")

        color_mapping = {(0.10,0.25):'blue',(0.25,0.50):'green',(0.50,0.75):'red',(0.75,1.00):'orange',(1.00,1.25):'purple',}

        print("Plotting all max amps and distance now........")

        # Create a figure and axis
        fig, ax = plt.subplots()
        # Loop through unique frequency bands
        for freq_range, group in dataframe.groupby('frequency band'):
            color = color_mapping[tuple(freq_range)]  # Get color based on frequency range
            ax.scatter(group['distance'], np.log10(group['max amplitude']), label = f'{freq_range[0]:.2f}-{freq_range[1]:.2f} Hz', color=color, marker = 'o',facecolors='none')
        # Add labels, legend, and show the plot
        plt.title('Original Data')
        ax.set_xlabel('Distance (km)')
        ax.set_ylabel('Max amplitude (Log scale)')
        plt.legend()
        #image = plt.savefig('AllMaxAmps.png') 
        image_filename = os.path.join('plots', 'AllMaxAmps.png')
        plt.savefig(image_filename)
        #save plot by date
        plt.show()  
    
    except: 
        print("Could not read Results.pkl or plot figure")