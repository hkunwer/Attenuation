import os
from locale import setlocale
import numpy as np
import csv
from numpy import abs
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

def check_existing_dataframe(): #checks if previous dataframe exists, or creates new dataframe if non-existent, This is to prevent overwriting data issues we had
    try:
        if os.path.exists('near_field_df.pkl'):
            with open('near_field_df.pkl', 'rb') as pkl_file:
                normalized_dataframe = pickle.load(pkl_file)
            print("Successful Check: File near_field_df.pkl exists")
            
        else:
            # If the file doesn't exist, create an empty dataframe
            print("Failed Check: Could not find existing dataframe")
            user_input = input(f"Do you want to save the data to a different file? (yes/no): ").lower()

            if user_input == 'yes':
                
                filename = input("Enter a new filename to save the data to (including extension, e.g., 'my_data.pkl'): ")
                                
                normalized_dataframe = pd.DataFrame(columns=["station name", "distance", "event", "location", 
                                                             "magnitude", "mag type", "max amplitude", "frequency band", 
                                                             "max amplitude (normalized)"])
                normalized_dataframe.to_pickle(filename)
                print("Completed: Created new dataframe to store normalized data later")
                
            else:
                print("Operation canceled.")
    
    except Exception as e:
        print(f"Error: {e}")
    
    return

def processANSStxt():

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
        
def filteringst(stream):
    
    taper=0.05
    stp = stream.copy()  # create backup
    # process data (maybe make this into a function)
    stp.detrend(type='polynomial', order=5) # pre-instrument removal
    stp.taper(taper)
    stp.remove_response(output="VEL") #using for comparison, will maybe use VEL for energies later. 
    stp.detrend(type='polynomial', order=5) # post-instrument removal
    stp.taper(taper)
 
    # Square the velocity values
    #stp.data = stp.data**2
    print("Completed: instrument response removal and taper, data returned as Velocity")
    
    return stp

def maxamp_calc_freq_bands(stream, EQ, Defaults, etime, eloc):
    
    #print("beginning max amp calculations")
    freqlist = [(0.10, 0.25),(0.25, 0.50),(0.50, 0.75),(0.75, 1.00),(1.00, 1.25)]
    
    stations = []
    distances = []
    event_dates = []
    event_location = []
    magnitudes = []
    magnitude_types = []
    max_amplitudes = []
    frequency_bands = []
    normalized_data = []
    
    for freq_range in freqlist:
        
        stp_freq = stream.copy()
        
        frequencymin, frequencymax = freq_range
        stp_freq.filter("bandpass", freqmin=frequencymin, freqmax=frequencymax)
        trace_count = 0  # Initialize a counter for the number of traces
        
        for tr in stp_freq:
            maxamp = np.max(abs(tr))
            threshold = 20 # Add a condition to exclude values above threshold
            if maxamp <= threshold:
                dist = tr.stats.distance / 1000
                station = tr.stats.station
                stations.append(station)
                distances.append(dist)
                event_dates.append(etime)
                event_location.append(eloc)
                magnitudes.append(EQ.Mag)
                magnitude_types.append(EQ.Mtype)
                max_amplitudes.append(maxamp)
                frequency_bands.append(freq_range)
                trace_count += 1 # Increment the trace counter
                
    print("Completed: frequency band filtering")
    
    if trace_count > 1:
        print("More than 1 trace detected, continuing")
        data = {
            'station name': stations,
            'distance': distances,
            'event': event_dates,
            'location':event_location,
            'magnitude': magnitudes,
            'mag type': magnitude_types,
            'max amplitude': max_amplitudes,
            'frequency band': frequency_bands,
        }

        # Create a new DataFrame with the trace data collected
        df = pd.DataFrame(data)
        outlier_removed_df = remove_outliers(df)
        normalized_dataframe = normalization(outlier_removed_df)
        return

    else:
        print("Only one or no traces available for normalization.")
        return 

def remove_outliers(df):
    
    # Calculate Z-scores for log-transformed max amplitude data
    # print("Removing outliers")
    # Calculate Z-scores for log-transformed max amplitude data
    z_scores = np.abs(stats.zscore(np.log(df['max amplitude'])))

    # Set a threshold for outlier detection (e.g., Z-score greater than 3)
    threshold = 4

    # Remove outliers based on the threshold
    outlier_removed_df = df[(z_scores < threshold)]

    return outlier_removed_df  

def normalization(outlier_removed_df):  # normalizes data by calling normalization factor calc and normalized max amp calc
    try:
        
        # Calculate the normalization factor (max amps) at a distance of 100 using the best-fit line for this event
        filtered_dataframe, normalization_factor = normalization_factor_calc(outlier_removed_df)
        # Apply the normalization function to create the 'max amplitude (normalized)' column
        outlier_removed_df['max amplitude (normalized)'] = outlier_removed_df.apply(
            lambda row: normalized_max_amplitude_calc(row, normalization_factor), axis=1
        ) 
        print("Completed: normalization")
        normalized_dataframe = outlier_removed_df.copy()
        
        # Save the updated dataframe to a pickle file
        save_dataframe_to_file(normalized_dataframe, filename="near_field_df.pkl")
        
        return normalized_dataframe

    except Exception as e:
        print(f"Incomplete: Could not complete normalization. Error: {e}")
        return outlier_removed_df
    
def normalization_factor_calc(df): #Calculates normalization factor for each event individually, using 1.0-1.25Hz freq band
    
    # Filter only the frequency band (1.0, 1.25) as a reference
    normalizeBench = 300
    filtered_dataframe = df[df['frequency band'] == (1.0, 1.25)]
    # Calculate the best fit line for the log-transformed filtered data
    filtered_log_max_amplitude = np.log10(filtered_dataframe['max amplitude'])
    filtered_coefficients = np.polyfit(filtered_dataframe['distance'], filtered_log_max_amplitude, 1)
    filtered_best_fit_line = np.poly1d(filtered_coefficients)
    
    # normalization factor
    normalization_factor = filtered_best_fit_line(normalizeBench)
    
    return filtered_dataframe, normalization_factor

def normalized_max_amplitude_calc(row, normalization_factor): #Calculates normalized max amps
    return np.log10(row['max amplitude']) - normalization_factor
        
def save_dataframe_to_file(dataframe, filename="near_field_df.pkl"):
    try:
        if os.path.exists(filename):
            # Load existing DataFrame from the file  
            existing_dataframe = pd.read_pickle(filename)
            print("Completed: Loading existing dataframe")
            # Append the new data to the existing DataFrame
            combined_dataframe = existing_dataframe.append(dataframe, ignore_index=True)

            # Save the combined DataFrame back to the file
            combined_dataframe.to_pickle(filename)

            print(f"Completed: Appended data to {filename}")
        else:
            dataframe.to_pickle(filename)
            print(f"Completed: Saved data to {filename}")
    except Exception as e:
        print(f"Error: {e}")
        
def event_normalization_plot(normalized_dataframe, event_name): #Create plots showing all traces and their distribution across different freq bands
    
#     #Create the "plots" folder if it doesn't exist
#     plots_folder = "/Users/hkunwer/Documents/research/EQenergy/Wiggles/Attenuation/plots"
#     print("testing")
#     if not os.path.exists(plots_folder):
#         os.mkdir(plots_folder)
#         print("Setting folder to plots")
    
    try: 
#         with open('near_field_df.pkl', 'rb') as pkl_file:
#             normalized_dataframe = pickle.load(pkl_file)
#             print("Completed: Loaded data from the near_field_df.pkl file") 

#         Drop rows where all columns are identical incase duplicates
#         normalized_dataframe = normalized_dataframe.drop_duplicates()
#         print("dropped duplicates")
        
        filtered_dataframe, filtered_best_fit_line = normalization_factor_calc(normalized_dataframe)    
        #print("filtered data")
        
        #Plotting
        #print("Plotting figures now")
        # Plot the log-transformed original data
        fig, axs = plt.subplots(1, 3, figsize=(18, 4))  
        # Subplot 1: Filtered data with best fit line
        normalizeBench = 100
        star_x = normalizeBench
        star_y = filtered_best_fit_line(star_x)
        filtered_log_max_amplitude = np.log10(filtered_dataframe['max amplitude'])
        
        axs[0].scatter(star_x, star_y, color='red', marker='*', s=50, label='Normalization Factor')
        axs[0].scatter(filtered_dataframe['distance'], filtered_log_max_amplitude, marker='o', label='(1.0,1.25)Hz', alpha=0.5)
        axs[0].plot(filtered_dataframe['distance'], filtered_best_fit_line(filtered_dataframe['distance']), color='red', label='Best Fit Line')
        axs[0].set_title('Filtered Data')
        axs[0].set_xlabel('Distance')
        axs[0].set_ylabel('Log10(Max Amplitude)')
        #print("Completed: plotting filtered dataframe")

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
            axs[1].scatter(group['distance'], np.log10(group['max amplitude']), marker='o', facecolors='none', alpha=0.5, label=f'{freq_range[0]:.2f}-{freq_range[1]:.2f} Hz', color=color)
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
            axs[2].scatter(group['distance'], group['max amplitude (normalized)'], marker='o', facecolors='none', alpha=0.5, label=f'{freq_range[0]:.2f}-{freq_range[1]:.2f} Hz', color=color)
        axs[2].set_title('Normalized Data')
        axs[2].set_xlabel('Distance')
        axs[2].set_ylabel('Log10(Max Amplitude) (Normalized)')
        axs[2].legend()
        
        # Set the plot title based on the 'event_name'
        plt.suptitle(f'Event: {event_name}')

        # Save the plot with a dynamic filename
        image_filename = os.path.join('plots', f'{event_name}_maxampPlot.png')
        print(f"Completed: saved {event_name}_maxampPlot.png")
        plt.savefig(image_filename)
        plt.tight_layout()
        plt.show()

    except: 
        print("Incomplete: Could not read near_field_df.pkl or plot figure")
        
def max_and_normalized_max_plots(df): #Plots original max amps against distance and compares to normalized max amps vs distance.

    freqlist = [(0.10, 0.25),(0.25, 0.50),(0.50, 0.75),(0.75, 1.00),(1.00, 1.25)]

    for freq_band in freqlist:

        filtered_df = df[df['frequency band'] == freq_band]
        #print(filtered_df)

        # Convert 'event' column to string type
        filtered_df.loc[:, 'event'] = filtered_df['event'].astype(str)

        # Get unique events (dates) from the filtered data
        unique_events = filtered_df['event'].unique()
        #print(unique_events)
        # Create a figure with two subplots
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))
        #print("plotting now")
        for event in unique_events:
            event_data = filtered_df[filtered_df['event'] == event]
            axs[0].scatter(
                event_data['distance'],
                np.log10(event_data['max amplitude']),
                label=event)
            axs[0].set_xlabel('Distance')
            axs[0].set_ylabel('Max amplitude')
            axs[0].set_title(f'Max amplitude for Frequency Band {freq_band}')
            axs[0].legend()
            axs[0].grid(True)

            axs[1].scatter(event_data['distance'],
                event_data['max amplitude (normalized)'],
                label=event)
            axs[1].set_xlabel('Distance')
            axs[1].set_ylabel('Max amplitude (normalized)')
            axs[1].set_title(f'Max amplitude (normalized) for Frequency Band {freq_band}')
            axs[1].legend()
            axs[1].grid(True)

        plt.tight_layout()
        plt.show()