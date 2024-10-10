# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 12:52:47 2018
Lecture 7.3/9
Reads several FIT-file from directory and plots 2D spectrum with background subtracted
Left mouse click saves x-/y-coords of a FIT-file 2D-image
Left mouse event.buttion  == 1
Right mouse event.button  == 3
Wheel/thumb event.buttion == 2
Click at least two points on the spectrum
@author: Christian Monstein, ETH Zurich
"""
from matplotlib import cm
import numpy as np
from matplotlib import pyplot as plt 
from astropy.io import fits
import glob
#-----------------------------------------------------------------------------------------------
# user input here:
#MyPath = 'LWA data/recent_4_23/' # Edit in case the file is elsewhere
#MyPath = 'LWA data/4_17/' # Edit in case the file is elsewhere
#MyPath = 'LWA data/4_19/' # Edit in case the file is elsewhere
MyPath = 'LWA data/4_22/'
#MyPath = 'LWA data/5_29/'
#MyPath = 'LWA data/5_21/'
#MyPath = 'LWA data/6_21/'

#MyFile = 'Beelink1_20240423_134500_59.fit'
#MyFile = 'ALASKA-COHOE_20240423_174459_62.fit' 
#MyFile = 'Beelink1_20240417_160000_59.fit' # file identification, use wildcards (*, ?)
#MyFile = 'Beelink1_20240419_134500_59.fit' # file identification, use wildcards (*, ?)
MyFile = 'Beelink1_20240422_171500_59.fit' 
#MyFile = 'UMICH_20240529_104500_59.fit'
#MyFile = 'UMICH_20240621_194500_59.fit'
# MyFile = 'UMICH_20240521_170000_59.fit'
harmonic = 1 # select excitation mode: 1=fundamental, 2=1st harmonic
newkirk  = 1.8 # select Newkirk model: 1.....4
zoomparameter = [0,900,10,85] # Starttime [s], Stoptime [s], Lowfrequency [MHz], Highfrequency [MHz]y
vmin = -5 # -5.....1 minimum visible value
vmax =  30 # 4....60 maximum visible value
# play with above parameters to get best quality and zoom of the burst
#-----------------------------------------------------------------------------------------------
Rsun = 695700. # Sun radius [km]

#--------------------------------------------------------------------------------------
paths = glob.glob(MyPath+MyFile)
print ('\nFIT-files found: ',(len(paths)))
#print(paths)

#-------------------------------------------------------------------------------
def get_data_from_fits(path):
    with fits.open(path) as hdu:
        data = hdu[0].data.astype(np.float32) # float16
    return data

def standardize_rows(arr):
    '''Standardize each row (subtract mean and divide by stdevs).'''
    mean_per_row = np.mean(arr, axis=1, keepdims=True)
    std_per_row = np.std(arr, axis=1, keepdims=True)
    epsilon = 1e-10
    std_per_row = np.maximum(std_per_row, epsilon)
    standardized_data = (arr - mean_per_row) / std_per_row
    return standardized_data

def normalize(arr):
    '''Normalize array to have values in range [0,255] for creating image.'''
    y = (arr - np.min(arr)) / (np.max(arr) - np.min(arr)) * 255
    return y


def remove_vertlines(data, threshold_low, threshold_high):
    """Remove spikes in the data by identifying and excluding time indices with spikes."""
    # denoised_data = []
    # median_values = np.median(data, axis=0)
    # spike_indices = np.where(median_values > threshold)[0]
    # for index in spike_indices:
    #     data[:, index] = np.median(data[:,0:10], axis=1)
    standardized = normalize(standardize_rows( data )).astype(np.uint8)
    median_values = np.median(standardized[:, :10], axis=1)     # Median value of the first 10 time columns, considered as background
    # Create a mask where data exceeds the upper threshold or is below the lower threshold. THIS IS GOOD AT REMOVING VERTICAL LINES, BUT IT ALSO ENLARGES OTHER NOISES
    mask = (standardized > threshold_high) | (standardized < threshold_low)
    # Replace values where the mask is True with the median values
    denoised_data = np.where(mask, median_values[:, np.newaxis], standardized)
    return denoised_data
#-------------------------------------------------------------------------------

instrument = fits.open(paths[0])[0].header['INSTRUME']
date       = fits.open(paths[0])[0].header['DATE-OBS']
# mytitle    = fits.open(paths[0])[0].header['CONTENT']
mytitle    = '2024/04/23  Radio flux density, e-CALLISTO (PM)'
mytitle    = '2024/04/23  Radio flux density, e-CALLISTO (COHOE)'
T0         = fits.open(paths[0])[0].header['TIME-OBS'] # start time as string
time       = fits.open(paths[0])[1].data[0][0] # time axis array
frequency  = fits.open(paths[0])[1].data[0][1] # frequency axis array

hh = float(T0.split(":")[0])
mm = float(T0.split(":")[1])
ss = float(T0.split(":")[2])
start_time = hh*3600 + mm*60 + ss # start time as float

#-------------------------------------------------------------------------------

data = np.hstack([get_data_from_fits(path) for path in paths])
print (data.shape)
dT = fits.open(paths[0])[0].header['CDELT1']
time_axis = (start_time + dT * np.arange(data.shape[1]))/3600.0
time_sec = (dT * np.arange(data.shape[1]))
dmean = np.mean(data, axis=1,keepdims=True)
dflat = data - dmean # subtraction background
#denoised_data = remove_vertlines(dflat, threshold_low=60, threshold_high=110) # For the type VI
denoised_data = remove_vertlines(dflat, threshold_low=50, threshold_high=130) # For the type II
#denoised_data = remove_vertlines(dflat, threshold_low=50, threshold_high=120) # For the type III
# denoised_data = remove_vertlines(dflat, threshold_low=60, threshold_high=160) # For the type V
# Subtract mean after vertical lines removement
dmean_2 = np.mean(denoised_data , axis=1,keepdims=True)
#denoised_data = denoised_data-dmean_2
# denoised_data = remove_vertlines(data, threshold_low=50, threshold_high=100)
#denoised_data = denoised_data-dmean
#-------------------------------------------------------------------------------

extent = (time_sec[0],time_sec[-1], frequency[-1],frequency[0]) # range for full 2D-image

fig, ax = plt.subplots(1,figsize=(9,7))
# im = ax.imshow(denoised_data, aspect = 'auto', extent = extent,cmap=cm.jet) # jet, CMRmap, gnuplot2, magma, plasma
# im = ax.imshow(dflat, aspect = 'auto', extent = extent,cmap=cm.jet) # jet, CMRmap, gnuplot2, magma, plasma
im = ax.imshow(dflat, aspect = 'auto', extent = extent,cmap=cm.CMRmap, norm=plt.Normalize(vmin, vmax)) # jet, CMRmap, gnuplot2, magma, plasma
cbar = plt.colorbar(im, ax=ax, label='Intensity [digit]')
# Adjust the fontsize of the colorbar
cbar.ax.tick_params(labelsize=35)  # Change the font size here
cbar.set_label('Intensity [digit]', fontsize=35)
T0 = str(13 + 4) + ":" + "45"
plt.xlabel("Time [s] after " + T0 + ' UT', fontsize=35)
plt.ylabel("Frequency [MHz]", fontsize=35)
plt.title(mytitle, fontsize=50)
plt.axis(zoomparameter)
plt.tick_params(axis='both', which='major', labelsize=30)
ax.autoscale = False # preventing plot from rescaling imag
#-----------------------------------------------------------------------------------------------

class MouseMonitor:
    fig = None
    axes = None

    def __init__(self, fig, ax):
        self.axes = ax
        self.fig = fig
        global coords, time, freq, Ne, rs, vr, dfdt
        coords = []
        time   = []
        freq   = []
        Ne     = []
        rs     = []
        vr     = []
        dfdt   = []         # This does not take intensity into account

        
    def __call__(self, event):
        if event.button == 3: # right mouse click terminates action
            fig.canvas.mpl_disconnect(cid)
            #plt.savefig(paths[0]+'.png')

            for i in range(len(coords)): # save all entries
                xn = coords[i][0] 
                yn = coords[i][1]
                time.append(xn)
                freq.append(yn)
                    
            for i in range(0,len(freq)):
                ne = (freq[i] / (harmonic*8.977e-3))**2.0 # electron density
                Ne.append(ne) # electron density
                rs.append(4.32 / (np.log10(ne/(newkirk*4.2e4)))) # radius scale
            
            plt.figure(figsize=(10,5))
            if (harmonic<2):
                st = ' , Newkirk model={:1.0f},'.format(newkirk) + ' fundamental'
            else:
                st = ' , Newkirk model={:1.0f},'.format(newkirk) + ' 1st harmonic'
            #plt.suptitle('CME speed analysis of ' + paths[0] + st, fontsize=30)
            tick_size = 20
            plt.suptitle('CME speed analysis of ' + '4/23 type II SRB' + st, fontsize=50)
            plt.subplot(2,3,1)
            plt.plot(time,freq,'-o',color="red")
            plt.grid()
            plt.xlabel("Time [s]",fontsize=25)
            plt.ylabel("Plasma frequency [MHz]",fontsize=25)
            plt.yticks(fontsize=tick_size)
            plt.xticks(fontsize=tick_size)
            
            plt.subplot(2,3,2)
            dfdt = np.abs(np.diff(freq) / np.diff(time))
            plt.plot(time[:-1],dfdt,'-o',color="green")
            print("Drift rate average:", np.mean(dfdt))         # Print out the average drift rate
            print("Drift rates:", dfdt) 
            plt.grid()
            plt.xlabel("Time [s]",fontsize=25)
            plt.ylabel("df/dt [MHz/s]",fontsize=25)
            plt.yticks(fontsize=tick_size)
            plt.xticks(fontsize=tick_size)
            
            plt.subplot(2,3,3)
            plt.plot(time,rs,'-o',color="blue")
            plt.grid()
            plt.xlabel("Time [s]",fontsize=25)
            plt.ylabel("Height [Rsun]",fontsize=25)
            plt.yticks(fontsize=tick_size)
            plt.xticks(fontsize=tick_size)
            
            plt.subplot(2,3,4)
            plt.plot(rs,freq,'-o',color="magenta")
            plt.grid()
            plt.ylabel("Plasma frequency [MHz]",fontsize=25)
            plt.xlabel("Height [Rsun]",fontsize=25)
            plt.yticks(fontsize=tick_size)
            plt.xticks(fontsize=tick_size)
    
            plt.subplot(2,3,5)
            plt.plot(rs[:-1],dfdt,'-o',color="cyan")
            plt.grid()
            plt.ylabel("Drift [MHz/s]",fontsize=25)
            plt.xlabel("Height [Rsun]",fontsize=25)
            plt.yticks(fontsize=tick_size)
            plt.xticks(fontsize=tick_size)
            
            plt.subplot(2,3,6)
            vr = np.diff(rs) / np.diff(time) * Rsun    # Velocity
            plt.plot(time[:-1],vr,'-o',color="black")
            plt.grid()
            plt.xlabel("Time [s]",fontsize=25)
            plt.ylabel("Speed [km/s]",fontsize=25)
            plt.yticks(fontsize=tick_size)
            plt.xticks(fontsize=tick_size)
            
            with open(paths[0]+'.table.txt', "w") as fp:   # Save x/y-data in file
                fp.write('    T[s],  F[MHz],   Ne[cm^-3], Rs[Rsun]\n') # write header information
                print ('    T[s],  F[MHz],   Ne[cm^-3], Rs[Rsun]')
                for i in range(0,len(freq)):
                    st = '{:8.3f},'.format(time[i]) + \
                         '{:8.3f},'.format(freq[i]) + \
                        '{:12.1f},'.format(Ne[i])   + \
                          '{:5.2f}'.format(rs[i])
                    fp.write(st+'\n')      
                    print (st)

            print ('\nStatistical results for CME velocity:')
            ym = np.mean(vr)        # Median value of all the radial velocity
            st = 'Mean      ={:7.1f}'.format(ym) + ' km/s'
            print (st)
            
            ym = np.median(vr)
            st = 'Median    ={:7.1f}'.format(ym) + ' km/s'
            print (st)
    
            v1 = (rs[-1] - rs[0]) / (time[-1] - time[0]) * Rsun  # rs: solar radius multiple
            st = '1st order ={:7.1f}'.format(v1) + ' km/s'
            print (st)
            #plt.savefig(paths[0]+'.results.png')


        else:    
            self.x = event.xdata
            self.y = event.ydata
            self.axes.plot(self.x, self.y, 'wo', linewidth = 10) #This don't work
            self.fig.canvas.draw_idle()
            #global coords
            coords.append((self.x, self.y)) # update x/y-array


mouse = MouseMonitor(fig, ax)

cid = fig.canvas.mpl_connect('button_press_event', mouse) 
#-----------------------------------------------------------------------------------------------
