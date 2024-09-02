from matplotlib import cm
import numpy as np
from matplotlib import pyplot as plt 
from astropy.io import fits
import glob

# User input here:
MyPath = 'LWA data/7_25/'
MyFile1 = 'ALASKA-COHOE_20240725_153000_63.fit'
MyFile2 = 'ALASKA-COHOE_20240725_154500_63.fit'


harmonic = 2  # Select excitation mode: 1=fundamental, 2=1st harmonic
newkirk  = 1.8  # Select Newkirk model: 1.....4
zoomparameter = [400,1200,30,85]  # Starttime [s], Stoptime [s], Lowfrequency [MHz], Highfrequency [MHz]y
vmin = -5  # Minimum visible value
vmax = 20  # Maximum visible value
Rsun = 695700.  # Sun radius [km]

# Get data from FITS files:
def get_data_from_fits(path):
    with fits.open(path) as hdu:
        data = hdu[0].data.astype(np.float32)  # float16
        print(f"Data shape from {path}: {data.shape}")
        print(f"Data min: {np.min(data)}, max: {np.max(data)}")
        print(f"Any NaNs? {np.isnan(data).any()}")
        print(f"Any zeros? {np.all(data == 0)}")
    return data

# Load both FITS files:
paths = [MyPath + MyFile1, MyPath + MyFile2]
data1 = get_data_from_fits(paths[0])
data2 = get_data_from_fits(paths[1])

# Check if data shapes match along the frequency axis
if data1.shape[0] != data2.shape[0]:
    raise ValueError("Frequency axes do not match between the two FITS files.")

# Concatenate data horizontally
data = np.concatenate((data1, data2), axis=1)

# Extract header information and time/frequency axes from both files:
def extract_header_info(path):
    with fits.open(path) as hdu:
        instrument = hdu[0].header['INSTRUME']
        date = hdu[0].header['DATE-OBS']
        mytitle = hdu[0].header['CONTENT']
        T0 = hdu[0].header['TIME-OBS']  # Start time as string
        time = hdu[1].data[0][0]  # Time axis array
        frequency = hdu[1].data[0][1]  # Frequency axis array
        dT = hdu[0].header['CDELT1']  # Time step
    return instrument, date, mytitle, T0, time, frequency, dT

# Extract info from both files
instrument1, date1, mytitle1, T01, time1, frequency1, dT1 = extract_header_info(paths[0])
instrument2, date2, mytitle2, T02, time2, frequency2, dT2 = extract_header_info(paths[1])

# Adjust the time axis for the second file
hh1, mm1, ss1 = map(float, T01.split(":"))
start_time1 = hh1 * 3600 + mm1 * 60 + ss1

hh2, mm2, ss2 = map(float, T02.split(":"))
start_time2 = hh2 * 3600 + mm2 * 60 + ss2

time_axis1 = (start_time1 + dT1 * np.arange(data1.shape[1]))/3600.0
time_axis2 = (start_time2 + dT2 * np.arange(data2.shape[1]))/3600.0
# Concatenate time axes
dT = fits.open(paths[0])[0].header['CDELT1']
time_axis = np.concatenate((time_axis1, time_axis2))
time_sec = (dT * np.arange(data.shape[1]))
# Background subtraction
dmean = np.mean(data, axis=1, keepdims=True)
dflat = data - dmean  # Subtract background
# Raw Data Plot
extent = (time_sec[0], time_sec[-1], frequency1[-1], frequency1[0])  # Range for full 2D-image
print(extent)
fig, ax = plt.subplots(1, figsize=(12, 7))
ax.imshow(dflat, aspect='auto', extent=extent, cmap=cm.CMRmap, norm=plt.Normalize(vmin, vmax))
plt.xlabel("Time [s] after " + T01 + ' UT', fontsize=14)
plt.ylabel("Frequency [MHz]", fontsize=14)
plt.title(mytitle1, fontsize=14)
plt.axis(zoomparameter)
plt.tick_params(axis='both', which='major', labelsize=14)
ax.autoscale = False  # Preventing plot from rescaling image

plt.show()

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
