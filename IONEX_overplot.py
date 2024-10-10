import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re
import cartopy.crs as ccrs
import wget
import subprocess
import matplotlib.cm
import matplotlib.colors
from datetime import datetime, timedelta
import matplotlib.dates as mdates

def parse_map(tecmap, exponent=-1):
    tecmap = re.split('.*END OF TEC MAP', tecmap)[0]
    return np.stack([np.fromstring(l, sep=' ') for l in re.split('.*LAT/LON1/LON2/DLON/H\\n', tecmap)[1:]]) * 10 ** exponent

def get_tecmaps(filename):
    with open(filename) as f:
        ionex = f.read()
        return [parse_map(t) for t in ionex.split('START OF TEC MAP')[1:]]

def get_tec(tecmap, lat, lon):
    i = round((87.5 - lat) * (tecmap.shape[0] - 1) / (2 * 87.5))
    j = round((180 + lon) * (tecmap.shape[1] - 1) / 360)
    return tecmap[i, j]

def ionex_filename(year, day, centre, zipped=True):
    return '{}g{:03d}0.{:02d}i{}'.format(centre, day, year % 100, '.Z' if zipped else '')

def ionex_ftp_path(year, day, centre):
    return 'ftp://cddis.gsfc.nasa.gov/gnss/products/ionex/{:04d}/{:03d}/{}'.format(year, day, ionex_filename(year, day, centre))

def ionex_local_path(year, day, centre='esa', directory='/tmp', zipped=False):
    return directory + '/' + ionex_filename(year, day, centre, zipped)

def download_ionex(year, day, centre='esa', output_dir='/tmp'):
    wget.download(ionex_ftp_path(year, day, centre), output_dir)
    subprocess.call(['gzip', '-d', ionex_local_path(year, day, centre, output_dir, zipped=True)])

def plot_tec_map(tecmap):
    proj = ccrs.PlateCarree()
    f, ax = plt.subplots(1, 1, subplot_kw=dict(projection=proj))
    ax.coastlines()
    h = plt.imshow(tecmap, cmap='viridis', vmin=0, vmax=100, extent=(-180, 180, -87.5, 87.5), transform=proj)
    plt.title('VTEC map')
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size='5%', pad=0.1, axes_class=plt.Axes)
    f.add_axes(ax_cb)
    cb = plt.colorbar(h, cax=ax_cb)
    plt.rc('text', usetex=True)
    cb.set_label('TECU ($10^{16} \\mathrm{el}/\\mathrm{m}^2$)')

day_cmap = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=1, vmax=365), cmap='inferno_r')

# Define the folder containing the INX files
folder_path = "Code/IONEX files/Eclipse_adjacent"

# Initialize an empty list to store TEC maps from all files
all_tecmaps = []

# Read and parse each INX file in the folder
for filename in os.listdir(folder_path):
    if filename.endswith(".INX"):
        inx_file = os.path.join(folder_path, filename)
        tecmaps = get_tecmaps(inx_file)
        all_tecmaps.append(tecmaps)  # Store each day's tecmaps separately

# Initialize an empty list to store TEC values
tec_values = []

# Get TEC values for each time interval
for day_tecmaps in all_tecmaps:
    day_tec_values = []
    for tecmap in day_tecmaps:
        day_tec_values.append(get_tec(tecmap, 40.46, -85.52))
    tec_values.append(day_tec_values)

# Create the x-axis values
time_intervals = np.arange(0, len(tec_values[0]) * 0.25, 0.25)

# Plot the combined TEC maps with different colors for each day
colors = ['b', 'g', 'r', 'c', 'm']
start_time = datetime(2024, 4, 6, 0, 0)
for day, day_tec_values in enumerate(tec_values):
    plt.plot(time_intervals, day_tec_values, color=colors[day % len(colors)], alpha=0.7, label=f'{(start_time + timedelta(days=day)).strftime("%Y/%m/%d")}')
    # if day < len(tec_values) - 1:
    #     # Connect the last point of the current day to the first point of the next day
    #     plt.plot([time_intervals[-1], time_intervals[0]], 
    #              [day_tec_values[-1], tec_values[day + 1][0]], 
    #              color=colors[(day + 1) % len(colors)], alpha=0.7)

plt.title('TEC variation of five consecutive days near 4/8 (Upland, Indiana) every 15mins')  # Indiana time zone = UTC-4
plt.xlabel('UTC hour')
plt.ylabel('VTEC (TECU)')

# Define the specific times for annotations
eclipse_start = 13.87 + 4
total_eclipse = 15.15 + 4
eclipse_end =   16.43 + 4

# Add annotations for eclipse events
plt.axvline(x=eclipse_start, color='red', linestyle='--')
plt.axvline(x=total_eclipse, color='red', linestyle='--')
plt.axvline(x=eclipse_end, color='red', linestyle='--')

plt.text(eclipse_start, 25, 'Eclipse start', fontsize=8, color='red', verticalalignment='top', horizontalalignment='center')
plt.text(total_eclipse, 24, 'Total eclipse', fontsize=8, color='red', verticalalignment='top', horizontalalignment='center')
plt.text(eclipse_end, 23, 'Eclipse end', fontsize=8, color='red', verticalalignment='top', horizontalalignment='center')

plt.grid(True)

# Format the x-axis to show UTC hours
# plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=4))
# plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# plt.gcf().autofmt_xdate()

# Add x-axis ticks every hour and label every 5 hours
hours = np.arange(0, 25, 1)
hour_labels = [f'{int(h)}:00' if h % 5 == 0 else '' for h in hours]
plt.xticks(hours, hour_labels)

plt.legend(title="Date")
plt.show()
