from astropy.table import Table
from astropy.table import join
import pandas as pd
import os
import shutil
import glob

from pysyd import utils
from pysyd.target import Target
from pysyd import plots
from astropy.stats import mad_std

tics = Table.read("dr19_tics.fits", format='fits')

names = []
for x in list(glob.glob('../../../blue/jtayar/cmarasco/dr19/pysyd_outputs/*')):
    names.append(x[49:-4])

failed = pd.read_csv('../../../blue/jtayar/cmarasco/dr19/failed.txt', names=['name'])
failed_names = [name[4:] for name in failed['name']]

def save_star(name, star):
    
    e_numax = mad_std(star.params['results']['parameters']['numax_smooth'])
    e_dnu = mad_std(star.params['results']['parameters']['dnu'])
    
    with open('../../../blue/jtayar/cmarasco/dr19/data/'+str(name)+'_cor_PS.txt', 'w') as f:
        freq = star.params['plotting']['parameters']['zoom_freq']
        power = star.params['plotting']['parameters']['zoom_pow']
        for i in range(0,len(freq)):
            f.write(str(freq[i])+'\t'+str(power[i])+'\n')
    
    # Update results file
    file_path = '../../../blue/jtayar/cmarasco/dr19/pysyd_outputs/'+str(name)+'.txt'
    tic_id = 'TIC ' + str(name)
    
    new_data = pd.DataFrame([[tic_id, star.params['numax_smoo'], e_numax, star.params['obs_dnu'], e_dnu]],
                            columns=['TIC ID', 'numax', 'e_numax', 'dnu', 'e_dnu'])
    new_data.to_csv(file_path, sep='\t', index=False, header=False)

def run_pysyd(name,numax,dnu):
    
    # if os.path.isfile('../../../blue/jtayar/cmarasco/dr19/info/star_info.csv'):
    #     os.remove('../../../blue/jtayar/cmarasco/dr19/info/star_info.csv')
    # f = open('../../../blue/jtayar/cmarasco/dr19/info/star_info.csv', 'w')
    # f.write('star,seed')
    # f.close()
    
    if os.path.isfile('info/star_info.csv'):
        os.remove('info/star_info.csv')
    f = open('info/star_info.csv', 'w')
    f.write('star,seed')
    f.close()

    print('setting up')
    params = utils.Parameters()
    params.add_targets(stars=name)

    star = Target(name, params)

    star.params['numax'] = numax
    star.params['dnu'] = dnu

    star.params['verbose'] = True
    star.params['mc_iter'] = 50
    star.params['smooth_ps'] = 0
    star.params['show'] = False
    star.params['overwrite'] = True

    # star.params['lower_ps'] = numax*(1-0.5)
    # star.params['upper_ps'] = numax*(1+0.8)

    star.params['lower_ex'] = 1
    star.params['upper_ex'] = 400
    
    star.params['inpdir'] = '../../../blue/jtayar/cmarasco/dr19/data'
    star.params['infdir'] = '../../../blue/jtayar/cmarasco/dr19/info'
    star.params['outdir'] = '../../../blue/jtayar/cmarasco/dr19/results'

    print('loading...')
    star.load_data()
    try:
        print('processing')
        star.process_star()
            
    except Exception as e:
        print('*********************************')
        print(name)
        print(e)
        print('*********************************')

    try:
        save_star(name,star)
    except Exception as e:
        tic_id = 'TIC ' + str(name)

        print(tic_id + ' save failed')
        print(e)

        # Update failed file
        file_path = '../../../blue/jtayar/cmarasco/dr19/failed.txt'
        with open(file_path, 'a') as file:
            file.write(tic_id+'\n')
            file.close()
    
    results_folder='../../../blue/jtayar/cmarasco/dr19/results/'+str(name)+'/'
    if os.path.exists(results_folder):
        shutil.rmtree(results_folder)
    if os.path.exists('results/'+str(name)+'/'):
        shutil.move('results/'+str(name)+'/', results_folder)

for i,tic in enumerate(tics['tic']):
    if os.path.isfile('../../../blue/jtayar/cmarasco/dr19/data/'+str(tic)+'_LC.txt') and str(tic) not in names and str(tic) not in failed_names:
        print(i, ' found')
        print(tic,tics['numax'][i],tics['est_dnu'][i])
        run_pysyd(tic,tics['numax'][i],tics['est_dnu'][i])