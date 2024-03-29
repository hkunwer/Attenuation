import os
from locale import setlocale
import numpy as np
import csv
from numpy import abs
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from datetime import datetime

def processANSStxt():
    """
    Processes data from 'ANSS_data.txt', formats and saves the processed data to a CSV file named 'ANSS_processed_data.csv'.
    The function then reads this CSV file into a pandas DataFrame and returns it.

    Returns:
    - ANSS: DataFrame containing the processed data.
    """
    try:
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
        print("Successfully processed ANSS data")
        return ANSS
    
    except Exception as e:
        print(f"Error processing ANSS data: {e}")
        return None

def check_existing_dataframe():
    """
    Checks if a dataframe stored as 'near_field_df.pkl' exists. If it exists, loads it; 
    otherwise, prompts the user to create a new dataframe to avoid overwriting issues.
    
    Returns:
    - normalized_dataframe: The loaded or newly created dataframe, or None if operation is canceled or fails.
    """
    try:
        if os.path.exists('near_field_df.pkl'):
            with open('near_field_df.pkl', 'rb') as pkl_file:
                normalized_dataframe = pickle.load(pkl_file)
            print("Successful Check: File near_field_df.pkl exists")
        else:
            print("Failed Check: Could not find existing dataframe.")
            user_input = input("Do you want to save the data to a different file? (yes/no): ").lower()

            if user_input == 'yes':
                filename = input("Enter a new filename to save the data to (including extension, e.g., 'my_data.pkl'): ")
                normalized_dataframe = pd.DataFrame(columns=["station name", "distance", "event", "location", 
                                                             "magnitude", "mag type", "frequency band", "max amplitude",
                                                             "max amplitude (normalized)"])
  
                normalized_dataframe.to_pickle(filename)
                print(f"Completed: Created new dataframe '{filename}' to store normalized data later.")
            else:
                print("Operation canceled.")
                return None  # Explicitly return None for clarity
    except Exception as e:
        print(f"Error encountered: {e}")
        return None  # Return None in case of error
    
    return normalized_dataframe    
    
def filter_stream_for_velocity_squared(stream):
    """
    Filters an ObsPy Stream to output velocity squared data.
    
    This function performs the following operations on the stream:
    1. Creates a backup copy of the original stream.
    2. Detrends the stream data using a 5th order polynomial, both before and after instrument response removal.
    3. Applies a taper to the stream with a specified percentage before and after instrument response removal.
    4. Removes the instrument response to output velocity data.
    5. Squares the velocity data to focus on energy calculations.
    
    Parameters:
    - stream: ObsPy Stream object to be processed.
    
    Returns:
    - Processed ObsPy Stream with velocity data squared.
    """
    
    taper_percentage = 0.05  # Taper percentage to use
    stp = stream.copy()  # Create a backup copy of the original stream
    
    # Pre-instrument response removal processing
    stp.detrend(type='polynomial', order=5)
    stp.taper(taper_percentage)
    
    # Remove instrument response, set to output velocity (VEL)
    stp.remove_response(output="VEL")
    
    # Post-instrument response removal processing
    stp.detrend(type='polynomial', order=5)
    stp.taper(taper_percentage)
    
    # Squaring the velocity data for energy calculation
    for trace in stp:
        trace.data = trace.data**2
    
    print("Completed: Instrument response removal, detrend, and taper. Data returned as Velocity Squared.")
    
    return stp

def maxamp_calc_freq_bands(stream, earthquake, defaults, event_time, event_location):
    """
    Calculates the maximum amplitude in specific frequency bands for a given seismic stream.
    
    Parameters:
    - stream: ObsPy Stream object containing seismic data.
    - earthquake: An object containing earthquake attributes like magnitude and type.
    - defaults: Not used in this snippet but could hold default configurations.
    - event_time: The time of the earthquake event.
    - event_location: The location of the earthquake event.
    
    Returns:
    - A DataFrame containing normalized data about max amplitudes in specified frequency bands.
    """
    
    freqlist = [(0.10, 0.25), (0.25, 0.50), (0.50, 0.75), (0.75, 1.00), (1.00, 1.25)]
    results = []

    for freq_range in freqlist:
        stp_freq = stream.copy()
        stp_freq.filter("bandpass", freqmin=freq_range[0], freqmax=freq_range[1])
        
        for tr in stp_freq:
            maxamp = np.max(np.abs(tr.data))
            #if maxamp <= 20:  # Assuming 20 is a meaningful threshold
            results.append({
                'station name': tr.stats.station,
                'distance': tr.stats.distance / 1000,
                'event': event_time,
                'location': event_location,
                'magnitude': earthquake.Mag,
                'mag type': earthquake.Mtype,
                'frequency band': freq_range,
                'max amplitude': maxamp,
 
            })

    df = pd.DataFrame(results)
    # Assuming remove_outliers and normalization are defined elsewhere
    outlier_removed_df = remove_outliers(df)
    normalization(outlier_removed_df)

def remove_outliers(df):
    """
    Removes outliers from a dataframe based on the Z-score of the log-transformed 'max amplitude' column.
    
    Parameters:
    - df: pandas DataFrame containing a 'max amplitude' column to be processed for outliers.
    
    Returns:
    - outlier_removed_df: DataFrame after outliers have been removed.
    """
    try:
        # Ensure 'max amplitude' column is in numeric format to avoid errors with log and Z-score calculations
        df['max amplitude'] = pd.to_numeric(df['max amplitude'], errors='coerce')

        # Calculate Z-scores for log-transformed 'max amplitude' data
        # Uses natural log; ensure no negative or zero values in 'max amplitude' to avoid NaN results
        z_scores = np.abs(stats.zscore(np.log(df['max amplitude'][df['max amplitude'] > 0])))

        # Define a threshold for outlier detection, e.g., Z-score greater than 4
        threshold = 4

        # Remove outliers based on the threshold
        # Note: Only rows with 'max amplitude' > 0 and within the Z-score threshold are retained
        outlier_removed_df = df[(z_scores < threshold) & (df['max amplitude'] > 0)]
        print("Completed: Outlier removal successful")
    except Exception as e:
        print(f"Error removing outliers: {e}")
        return df  # Return the original DataFrame in case of error

    return outlier_removed_df 

def normalization(outlier_removed_df):
    """
    Normalizes the maximum amplitude data using distance as the normalization factor.
    This involves calculating a normalization factor for the data based on a specified frequency band
    and applying this factor to each data point to create a normalized maximum amplitude column.

    Parameters:
    - outlier_removed_df: DataFrame containing outlier-removed data including 'distance' and 'frequency band' columns.

    Returns:
    - The DataFrame with an additional 'max amplitude (normalized)' column.
    """
    try:
        
        # Calculate the normalization factor (max amps) at a distance of 200 using the best-fit line for this event
        filtered_dataframe, normalization_factor = normalization_factor_calc(outlier_removed_df)
        
        # Apply the normalization function to create the 'max amplitude (normalized)' column
        outlier_removed_df['max amplitude (normalized)'] = outlier_removed_df.apply(
            lambda row: normalized_max_amplitude_calc(row, normalization_factor), axis=1
        ) 
        print("Completed: normalization")
        normalized_dataframe = outlier_removed_df.copy()
        
        # Save the updated dataframe to a pickle file
        save_dataframe_to_file(normalized_dataframe)

    except Exception as e:
        print(f"Incomplete: Could not complete normalization. Error: {e}")
    
def normalization_factor_calc(df):
    """
    Calculates the normalization factor for maximum amplitude based on distance.
    This is done by filtering the data for a specific frequency band and then calculating 
    the best-fit line for the log-transformed maximum amplitude data against v^2.

    Parameters:
    - df: DataFrame containing the data to normalize, including 'distance' and 'max amplitude' columns.

    Returns:
    - A tuple of the filtered DataFrame and the normalization factor calculated.
    """
    # Filter only the frequency band (1.0, 1.25) as a reference
    normalizeBench = 200
    filtered_dataframe = df[df['frequency band'] == (1.0, 1.25)]
    # Calculate the best fit line for the log-transformed filtered data
    filtered_log_max_amplitude = np.log10(filtered_dataframe['max amplitude'])
    filtered_coefficients = np.polyfit(filtered_dataframe['distance'], filtered_log_max_amplitude, 1)
    filtered_best_fit_line = np.poly1d(filtered_coefficients)
    
    # normalization factor
    normalization_factor = filtered_best_fit_line(normalizeBench)
    
    return filtered_dataframe, normalization_factor

def normalized_max_amplitude_calc(row, normalization_factor):
    """
    Calculates the normalized maximum amplitude for a single row/data point.

    Parameters:
    - row: A row from the DataFrame containing 'max amplitude' and 'distance' data.
    - normalization_factor: The calculated normalization factor to apply.

    Returns:
    - The normalized maximum amplitude value for the row.
    """
    return np.log10(row['max amplitude']) - normalization_factor

def save_dataframe_to_file(dataframe, directory="/Users/hkunwer/Documents/research/EQenergy/Attenuation/V2 Max Amps", base_filename="near_field_df.pkl"):
    """
    Saves a DataFrame to a file in the specified directory. If "near_field_df.pkl" exists in the directory, appends the new data to it.
    If the file doesn't exist, creates a new file with a timestamp to avoid potential overwriting.

    Parameters:
    - dataframe: The DataFrame to be saved or appended.
    - directory: The directory where the file is saved or to be saved. Defaults to the specified path.
    - base_filename: The base name of the file to save the DataFrame to. Defaults to 'near_field_df.pkl'.
    """
    if dataframe.empty:
        print("The DataFrame is empty. No data to save.")
        return

    # Generate the full file path
    file_path = os.path.join(directory, base_filename)

    try:
        # Check if the file already exists in the specified directory
        if os.path.exists(file_path):
            # Load the existing DataFrame and append the new data
            existing_dataframe = pd.read_pickle(file_path)
            combined_dataframe = pd.concat([existing_dataframe, dataframe], ignore_index=True)
            # Save the combined DataFrame back to the original file
            combined_dataframe.to_pickle(file_path)
            print(f"Data appended to {file_path} successfully.")
        else:
            # If the file does not exist, save the new DataFrame as the file
            dataframe.to_pickle(file_path)
            print(f"File {file_path} not found. Created new file with provided dataframe.")

    except Exception as e:
        # For any other exceptions, provide an error message
        print(f"An error occurred while trying to save the dataframe: {e}.")

# def event_normalization_plot(normalized_dataframe, event_name): #Create plots showing all traces and their distribution across different freq bands
    
# #     #Create the "plots" folder if it doesn't exist
# #     plots_folder = "/Users/hkunwer/Documents/research/EQenergy/Wiggles/Attenuation/plots"
# #     print("testing")
# #     if not os.path.exists(plots_folder):
# #         os.mkdir(plots_folder)
# #         print("Setting folder to plots")
    
#     try: 
# #         with open('near_field_df.pkl', 'rb') as pkl_file:
# #             normalized_dataframe = pickle.load(pkl_file)
# #             print("Completed: Loaded data from the near_field_df.pkl file") 

# #         Drop rows where all columns are identical incase duplicates
# #         normalized_dataframe = normalized_dataframe.drop_duplicates()
# #         print("dropped duplicates")
        
#         filtered_dataframe, filtered_best_fit_line = normalization_factor_calc(normalized_dataframe)    
#         #print("filtered data")
        
#         #Plotting
#         #print("Plotting figures now")
#         # Plot the log-transformed original data
#         fig, axs = plt.subplots(1, 3, figsize=(18, 4))  
#         # Subplot 1: Filtered data with best fit line
#         normalizeBench = 100
#         star_x = normalizeBench
#         star_y = filtered_best_fit_line(star_x)
#         filtered_log_max_amplitude = np.log10(filtered_dataframe['max amplitude'])
        
#         axs[0].scatter(star_x, star_y, color='red', marker='*', s=50, label='Normalization Factor')
#         axs[0].scatter(filtered_dataframe['distance'], filtered_log_max_amplitude, marker='o', label='(1.0,1.25)Hz', alpha=0.5)
#         axs[0].plot(filtered_dataframe['distance'], filtered_best_fit_line(filtered_dataframe['distance']), color='red', label='Best Fit Line')
#         axs[0].set_title('Filtered Data')
#         axs[0].set_xlabel('Distance')
#         axs[0].set_ylabel('Log10(Max Amplitude)')
#         #print("Completed: plotting filtered dataframe")

#         # Add equation of the line as a text annotation (adjust position)
#         equation = f'Line Equation: y = {filtered_best_fit_line[0]:.3f}x + {filtered_best_fit_line[1]:.3f}'
#         axs[0].text(0.30, 0.10, equation, transform=axs[0].transAxes, fontsize=9)  # Adjust position here
#         axs[0].legend()

#         # Subplot 2: Original data with different colors for each frequency band (excluding 1.0,1.25 Hz)
#         color_mapping = {
#             (0.10, 0.25): 'blue',
#             (0.25, 0.50): 'green',
#             (0.50, 0.75): 'red',
#             (0.75, 1.00): 'orange',
#             (1.00, 1.25): 'purple',
#         }
#         for freq_range, group in normalized_dataframe.groupby('frequency band'):
#             color = color_mapping[tuple(freq_range)]
#             axs[1].scatter(group['distance'], np.log10(group['max amplitude']), marker='o', facecolors='none', alpha=0.5, label=f'{freq_range[0]:.2f}-{freq_range[1]:.2f} Hz', color=color)
#         axs[1].set_title('Original Data')
#         axs[1].set_xlabel('Distance')
#         axs[1].set_ylabel('Log10(Max Amplitude)')
#         axs[1].legend()

#         # Subplot 3: Normalized data with different colors for each frequency band (including 1.0,1.25 Hz)

#         star_x = normalizeBench
#         star_y = filtered_best_fit_line(star_x) - filtered_best_fit_line(star_x)
#         axs[2].scatter(star_x, star_y, color='red', marker='*', s=50, label='Normalization Factor')

#         for freq_range, group in normalized_dataframe.groupby('frequency band'):
#             color = color_mapping[tuple(freq_range)]
#             axs[2].scatter(group['distance'], group['max amplitude (normalized)'], marker='o', facecolors='none', alpha=0.5, label=f'{freq_range[0]:.2f}-{freq_range[1]:.2f} Hz', color=color)
#         axs[2].set_title('Normalized Data')
#         axs[2].set_xlabel('Distance')
#         axs[2].set_ylabel('Log10(Max Amplitude) (Normalized)')
#         axs[2].legend()
        
#         # Set the plot title based on the 'event_name'
#         plt.suptitle(f'Event: {event_name}')

#         # Save the plot with a dynamic filename
#         image_filename = os.path.join('plots', f'{event_name}_maxampPlot.png')
#         print(f"Completed: saved {event_name}_maxampPlot.png")
#         plt.savefig(image_filename)
#         plt.tight_layout()
#         plt.show()

#     except: 
#         print("Incomplete: Could not read near_field_df.pkl or plot figure")
        
# def max_and_normalized_max_plots(df): #Plots original max amps against distance and compares to normalized max amps vs distance.

#     freqlist = [(0.10, 0.25),(0.25, 0.50),(0.50, 0.75),(0.75, 1.00),(1.00, 1.25)]

#     for freq_band in freqlist:

#         filtered_df = df[df['frequency band'] == freq_band]
#         #print(filtered_df)

#         # Convert 'event' column to string type
#         filtered_df.loc[:, 'event'] = filtered_df['event'].astype(str)

#         # Get unique events (dates) from the filtered data
#         unique_events = filtered_df['event'].unique()
#         #print(unique_events)
#         # Create a figure with two subplots
#         fig, axs = plt.subplots(1, 2, figsize=(12, 5))
#         #print("plotting now")
#         for event in unique_events:
#             event_data = filtered_df[filtered_df['event'] == event]
#             axs[0].scatter(
#                 event_data['distance'],
#                 np.log10(event_data['max amplitude']),
#                 label=event)
#             axs[0].set_xlabel('Distance')
#             axs[0].set_ylabel('Max amplitude')
#             axs[0].set_title(f'Max amplitude for Frequency Band {freq_band}')
#             axs[0].legend()
#             axs[0].grid(True)

#             axs[1].scatter(event_data['distance'],
#                 event_data['max amplitude (normalized)'],
#                 label=event)
#             axs[1].set_xlabel('Distance')
#             axs[1].set_ylabel('Max amplitude (normalized)')
#             axs[1].set_title(f'Max amplitude (normalized) for Frequency Band {freq_band}')
#             axs[1].legend()
#             axs[1].grid(True)

#         plt.tight_layout()
#         plt.show()