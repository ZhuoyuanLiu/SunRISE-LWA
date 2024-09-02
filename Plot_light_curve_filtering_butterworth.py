'''
    Time range: 3 hours 50 mins
    freq range: 45-870MHz
    
    Choose 45MHz, 100Hz, 160MHz to plot light curves after std subtraction 
 '''
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import math
from matplotlib import cm
# Path to the directory containing the FIT files
directory_path = 'Dataset\INDIANA'
#directory_path = 'Dataset\Peach Mountain'

# List to hold data arrays from all files
spectrogram_data = []
time_data = []
frequency_data = []


# Read each file in the directory
for filename in os.listdir(directory_path):
    if filename.endswith('.fits') or filename.endswith('.fit'):
        # Full path to the file
        filepath = os.path.join(directory_path, filename)
        
        # Open the FITS file
        with fits.open(filepath) as hdul:
            # Assuming the spectrogram data is in the first extension
            data = hdul[0].data
            header = hdul[0].header
    

        # Extract time and frequency information
            time_start_str = header['DATE-OBS'] + 'T' + header['TIME-OBS']
            time_start = datetime.strptime(time_start_str, '%Y/%m/%dT%H:%M:%S.%f')
            time_interval = header['CDELT1']  # Time interval per pixel (in seconds)
            num_time_points = data.shape[1]

            # Generate time data for this file
            if len(time_data) > 0:
                time_start = time_data[-1] + timedelta(seconds=time_interval)
            file_time_data = [time_start + timedelta(seconds=i*time_interval) for i in range(num_time_points)]

            if not spectrogram_data:
                # Extract frequency information only once (assuming it is the same for all files)
                # freq_start = header['CRVAL2']
                freq_start = 45
                # freq_interval = header['CDELT2']  # Frequency interval per pixel (in MHz)
                channel_num = header['CRVAL2']
                freq_interval = (90-45)/channel_num  # Frequency interval per pixel (in MHz)
                frequency_data = [freq_start + i*freq_interval for i in range(data.shape[0])]

            spectrogram_data.append(data)
            time_data.extend(file_time_data)
# Concatenate all data arrays along the time axis (assuming time is the second axis)
combined_spectrogram_orig = np.concatenate(spectrogram_data, axis=1)
# TO PLOT from 0 to 160MHz, can take the first valid channels

valid_channels = math.ceil((90-45)/freq_interval)
# valid_channels = -1
print(valid_channels)
combined_spectrogram = combined_spectrogram_orig[-valid_channels:-1]     # Slice from 45-160MHz, now the order in it is 160-45.
#combined_spectrogram = combined_spectrogram_orig 

# Plotting the combined spectrogram
plt.figure(figsize=(10, 5))
ax = plt.gca()
im = ax.imshow(combined_spectrogram, aspect='auto', origin='upper', extent=[mdates.date2num(time_data[0]), mdates.date2num(time_data[-1]), frequency_data[0], frequency_data[valid_channels-1]], cmap=cm.CMRmap)
ax.xaxis_date()
date_format = mdates.DateFormatter('%H:%M:%S')
ax.xaxis.set_major_formatter(date_format)
# Customize the time tick intervals
ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=15))  # Major ticks every minute
ax.xaxis.set_minor_locator(mdates.SecondLocator(interval=30))  # Minor ticks every 15 seconds

# Define the specific times for annotations
eclipse_start = datetime.strptime("2024/04/08T13:45:00", "%Y/%m/%dT%H:%M:%S")
total_eclipse = datetime.strptime("2024/04/08T15:10:00", "%Y/%m/%dT%H:%M:%S")
eclipse_end = datetime.strptime("2024/04/08T16:30:00", "%Y/%m/%dT%H:%M:%S")
# Add annotations for eclipse events
plt.axvline(x=mdates.date2num(eclipse_start), color='red', linestyle='--')
plt.axvline(x=mdates.date2num(total_eclipse), color='red', linestyle='--')
plt.axvline(x=mdates.date2num(eclipse_end), color='red', linestyle='--')

plt.text(mdates.date2num(eclipse_start), combined_spectrogram.shape[0], 'Eclipse start', color='red', verticalalignment='top', horizontalalignment='center')
plt.text(mdates.date2num(total_eclipse), combined_spectrogram.shape[0], 'Total eclipse', color='red', verticalalignment='top', horizontalalignment='center')
plt.text(mdates.date2num(eclipse_end), combined_spectrogram.shape[0], 'Eclipse end', color='red', verticalalignment='top', horizontalalignment='center')

plt.colorbar(im, label='Intensity')
plt.title('Combined Spectrogram of Original Data')
plt.xlabel('Time [HH:MM:SS]')
plt.ylabel('Frequency [MHz]')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()


# Plotting the standard deviations as a histogram
channel_std = np.std(combined_spectrogram, axis=1)
channel_mean = np.mean(combined_spectrogram, axis=1)
# plt.figure(figsize=(10, 5))
# plt.bar( frequency_data[:valid_channels], channel_std)  # Channel numbers starting from 1
# plt.title('Standard Deviation per Frequency Channel')
# plt.xlabel('Frequency Channel Value(MHz)')
# plt.ylabel('Standard Deviation of intensity')
# plt.grid(True)
# plt.show()


# Substract std from image
channel_std = np.repeat(channel_std[:, np.newaxis], combined_spectrogram.shape[1], axis=1)
channel_mean = np.repeat(channel_mean[:, np.newaxis], combined_spectrogram.shape[1], axis=1)
# sub_bg_img = combined_spectrogram-channel_std
sub_bg_img = combined_spectrogram-channel_mean
# Define the desired frequencies
desired_frequencies = [52.5, 55, 62.5, 67.5, 70, 82.5]       # Too many outliers at 45 and 160. For Indiana, looks like above 155MHz there is interesting sig, for PM, not easy to say
# desired_frequencies = [70, 74, 80, 86, 94, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 153 ]      # Stable set Ind. 74MHz was reported in E-call website with trend but no in mine
# desired_frequencies = [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800 ]  
# desired_frequencies = [205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285 ]  
# Find the indices of the desired frequencies
desired_indices = [np.argmin(np.abs(np.array(frequency_data) - freq)) for freq in desired_frequencies]      # Frequency_data is from 45-870MHz

#desired_indices = -desired_indices # 6/13: Update Order, reversed

print("desired_indices:", desired_indices)
# Ensure the time data is in numpy datetime64 format for consistent plotting
time_data_np = np.array( time_data, dtype='datetime64')
# Plotting the light curves

# Define a moving average filter function
def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')
window_size = 160       # Higher this value higher the smoothness

# Define a Savitzky-Golay smoothing filter
from scipy.signal import savgol_filter
import scipy.signal

cutoff = 0.0005 # Low-pass cutoff freq
sample_rate = int(1/time_interval)
sos = scipy.signal.butter(5, cutoff, 'lowpass', fs = sample_rate, output = 'sos')       # Construct butterworth filter.
                                                                                        # Fs is the sampling rate of the e-callisto signals, 4Hz
print( np.array(time_data).shape, sub_bg_img[1,:].shape)
plt.figure(figsize=(12, 6))

# Plot light curves
for i, freq in zip(desired_indices, desired_frequencies):
    # smoothed_curve = moving_average(sub_bg_img[i-1, :], window_size)
    smoothed_curve = scipy.signal.sosfiltfilt(sos, sub_bg_img[-i-1, :]) # Last one is 45, 
    #smoothed_curve = savgol_filter(sub_bg_img[i-1, :], 80, 2)
    smoothed_time = time_data_np[(window_size - 1) // 2 : -(window_size // 2)]
    # plt.plot( np.array(time_data), sub_bg_img[i-1, :], linestyle = '-', label=f'{freq} MHz')
    # This is the filter I used to smooth NASA eclipse X-band data
    #plt.plot( smoothed_time , smoothed_curve, linestyle = '-', label=f'{freq} MHz')
    plt.plot( time_data_np , smoothed_curve, linestyle = '-', label=f'{freq} MHz')

# Setting up date format and ticks
date_format = mdates.DateFormatter('%H:%M:%S')
plt.gca().xaxis.set_major_formatter(date_format)
plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(interval=15))  # Major ticks every 10 minutes
plt.gca().xaxis.set_minor_locator(mdates.MinuteLocator(interval=30))  # Minor ticks every 2 minutes

# Define the specific times for annotations
eclipse_start = datetime.strptime("2024/04/08T13:45:00", "%Y/%m/%dT%H:%M:%S")
total_eclipse = datetime.strptime("2024/04/08T15:10:00", "%Y/%m/%dT%H:%M:%S")
eclipse_end = datetime.strptime("2024/04/08T16:30:00", "%Y/%m/%dT%H:%M:%S")
# Add annotations for eclipse events
plt.axvline(x=mdates.date2num(eclipse_start), color='red', linestyle='--')
plt.axvline(x=mdates.date2num(total_eclipse), color='red', linestyle='--')
plt.axvline(x=mdates.date2num(eclipse_end), color='red', linestyle='--')

plt.text(mdates.date2num(eclipse_start), plt.ylim()[1] * 0.95, 'Eclipse start', color='blue', verticalalignment='top', horizontalalignment='center')
plt.text(mdates.date2num(total_eclipse), plt.ylim()[1] * 0.95, 'Total eclipse', color='blue', verticalalignment='top', horizontalalignment='center')
plt.text(mdates.date2num(eclipse_end), plt.ylim()[1] * 0.95, 'Eclipse end', color='blue', verticalalignment='top', horizontalalignment='center')

plt.xlabel('Time [HH:MM:SS]')
plt.ylabel('Intensity')
# plt.title(f'Using Butterworth filter with cutoff {cutoff} Hz, Smoothed Light Curves at Selected Frequencies (PM)')
plt.title(f'Using Butterworth filter with cutoff {cutoff} Hz, Smoothed Light Curves at Selected Frequencies (Indiana)')
plt.legend()
plt.grid(True)
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()