import os
import glob
import shutil
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from astropy.io import fits
from scipy.fft import fft, ifft
import pandas as pd
import astropy.units as u
import astropy.constants as c
from tqdm import tqdm

from scipy.fft import fft, fftfreq

import lightkurve as lk
from lightkurve import search_lightcurve


def lc_download(name, plot=True, overwrite=False, window=299, combine=True, min_sectors=0):
    if (not os.path.isfile('../../../blue/jtayar/cmarasco/dr19/data/'+name[4:]+'_LC.txt')) or overwrite or plot:
        # print('Downloading LC for:  '+name)
        # Search for light curves with the given name and author 'QLP'
        search_result = lk.search_lightcurve(name, author='QLP')
        # Extract the target name from the search result
        names = search_result.target_name
        for i in names:
            if i != names[0]:
                raise Exception("names in search don't match")
            else:
                name = names[0]

        # Download all light curves returned by the search
        lc_collect_og = search_result.download_all()
        # Create empty list for storing light curves
        lc_collect = []
        lc_collected = []
        lc_len = 0

        if len(lc_collect_og) >= min_sectors:
            if plot:
                fig, axs = plt.subplots(nrows=math.ceil(len(lc_collect_og)/4), ncols=4,
                                    figsize=[20, 1.5*math.ceil(len(lc_collect_og)/4)])
            for k in range(len(lc_collect_og)):
                if lc_collect_og[k].sector not in [8]:
                    lc_len += 1
                    lc_new = lc_collect_og[k][50:-50].normalize().flatten(window_length=window)
                    lc_new = lc_new.remove_outliers(sigma=2.75)
                    lc_collect.append(lc_new)

                    some_times = [time.value for time in lc_new['time']]
                    some_fluxes = [flux.value for flux in lc_new['flux']]

                    # Plot the light curve that will be saved
                    if plot:
                        if len(lc_collect_og) <= 4:
                            ax = axs[k]
                        else:
                            ax = axs[math.floor(k/4),k%4]
                        ax.plot(some_times,some_fluxes, color='C'+str(k))
                        # ax.legend(fontsize=28,borderpad=0,handlelength=0,loc='lower right')
                        txt = ax.text(min(some_times)+(0.85*(max(some_times)-min(some_times))),
                                min(some_fluxes)+(0.05*(max(some_fluxes)-min(some_fluxes))),
                                str(lc_collect_og[k].sector),
                                weight='black',size=32,color='k')
                        txt.set_path_effects([PathEffects.withStroke(linewidth=3,
                                            foreground='w')])


            lc_collected = lk.LightCurveCollection(lc_collect)
            lc = lc_collected.stitch().normalize()

            # Putting the LC data into text files that pysyd will be able to read
            times = lc["time"]
            fluxes = lc["flux"]
            sorted_indices = np.argsort(times.value)
            times = times[sorted_indices].value
            fluxes = fluxes[sorted_indices].value

            # Removing long gaps in time between sectors
            if combine:
                for i in range(1,len(times)):
                    if (times[i]-times[i-1]) > 10:
                        dif = times[i]-times[i-1]
                        times[i:] = times[i:]-dif

            if plot:
                plt.show()
                fig = plt.figure(figsize=[20, 1.5])
                plt.plot(times, fluxes, color='r')
                plt.show()
                
            # Saving the light curve to a txt file
            with open('../../../blue/jtayar/cmarasco/dr19/data/'+name+'_LC.txt', 'w') as f:
                for i in range(0,len(times)):
                    if "——" not in str(fluxes[i]):
                        f.write('\t'+str(times[i])+'\t'+str(fluxes[i])+'\n')
                f.close()

            # Updating the star_info table with new star
            with open('../../../blue/jtayar/cmarasco/dr19/info/star_info.csv', 'r') as f:
                if name not in f.read():
                    with open('../../../blue/jtayar/cmarasco/dr19/info/star_info.csv', 'a') as f:
                        f.write(name+',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'+ \
                                ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n')
                f.close()

            # Updating the sectors file with new star's number of sectors
            with open('../../../blue/jtayar/cmarasco/dr19/sectors.txt', 'r') as f:
                if name not in f.read():
                    with open('../../../blue/jtayar/cmarasco/dr19/sectors.txt', 'a') as f:
                        f.write(name+'\t'+str(lc_len)+'\n')
                f.close()

            # Remove power spectrum if previously calculated by pysyd
            if os.path.isfile('../../../blue/jtayar/cmarasco/dr19/data/'+name+'_PS.txt'):
                os.remove('../../../blue/jtayar/cmarasco/dr19/data/'+name+'_PS.txt')
                
                
with open('dr19_tics.txt', 'r') as f:
    names = f.readlines()
for i in range(len(names)):
    names[i] = 'TIC '+names[i].strip()
    
for name in tqdm(names):
    try:
        lc_download(name,plot=False)
    except Exception as e:
        print(e)
        
    # Clear cache
    cache = list(glob.glob('../.lightkurve/cache/mastDownload/HLSP/*'))
    for folder in cache:
        shutil.rmtree(folder)