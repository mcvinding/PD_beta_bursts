# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
from mne.preprocessing import create_ecg_epochs, create_eog_epochs, read_ica
from mne.io import read_raw_fif
from mne import pick_types
from os import listdir, chdir, mkdir, path
import os.path as op
import numpy as np
import csv
import sys
sys.path.append('/home/mikkel/PD_motor/global_scripts')
from PD_motor_pyfun import run_ICA_wMNE


# %% run options
project_part        = 'rest_ec'
overwrite_old_files = False     # Wheter files should be overwritten if already exist
do_run_ica          = True               #Whether tho run ICA (set to False if only testing)
sub_to_run          = 'all'          # A single string specifying a single subject to rerun or 'all'

#max95subs           = ['0319','0352','0314','0376','0377']   #Subjects with lowered Maxfilter correlation limit (i.e 0.95)

#FOLDERS ETC.
raw_path = '/home/mikkel/PD_motor/'+project_part+'/raw' # Files sorted in raw folder by means of symbolic links
meg_path = '/home/mikkel/PD_motor/'+project_part+'/meg_data' # Work directory for MEG
log_path = '/home/mikkel/PD_motor/'+project_part+'/summary'

# %% Run through files
allReadyRunFiles = []
file_list = listdir(raw_path)

if not sub_to_run=='all':
    if not any([f for f in file_list if any(xf in f for xf in sub_to_run)]):
#        print('No files for subject named: \"'+sub_to_run+'\"')
        file_list = []
    else:
        file_list = [f for f in file_list if any(xf in f for xf in sub_to_run)]
                    
## Find max95  subs            
#idx = [i for i,j in enumerate(file_list) if any(jx in j for jx in max95subs)]
#if idx:
#    file_list95 = [f for f in file_list if 'max0.95' in f]         #Use 95% mafiler cutoff data
#    file_list = [file_list[x] for x,y in enumerate(file_list) if x not in idx]
#    file_list = file_list+file_list95
#    file_list.sort()
    
#%%% Run
for ff in file_list:   
    sub = ff[:4]

    #Find all files for same subject to run ICA on.
    print('Now reading file: '+ff)
    prefx = ff[:12]
     
    if ff in allReadyRunFiles:
        print('...Nope! I already have done ICA on that! Moving to next file.')
        continue
    else:
        allReadyRunFiles += [f for f in file_list if f.startswith(prefx)]
        inFiles = [op.join(raw_path,f) for f in file_list if f.startswith(prefx)]
        print('Number of files to read = '+str(len(inFiles)))
        
    print('Now running ICA for sub = ' + ff[:4])

    # Make dirs, etc. for output
    sub_path = op.join(meg_path,sub)
    ica_path = op.join(meg_path,sub,'ica')
    out_icaFname = project_part+'-ica.fif'

    if not op.exists(sub_path):
        mkdir(sub_path)
    if not op.exists(ica_path):
        mkdir(ica_path)      
        
    if path.isfile(ica_path+'/'+out_icaFname):
        print('File "'+ica_path+'/'+out_icaFname+'" already exists')
        if not overwrite_old_files:
            print('Do not overwrite ICA file')
            continue
        else:
            print('Do overwrite ICA again! Running everything again!')
    
    #RUN ICA
    chdir(sub_path)
    if do_run_ica:
        run_ICA_wMNE(inFiles, ica_path, out_icaFname, prefilt=[1,40])  #,ica_fname=out_fname)
        print('--------------ICA DONE FOR '+sub+' ------------------') 
    else:
        print('-----------ICA *NOT* DONE FOR '+sub+' ---------------') 
            
#%% Find ECG and EOG artefacts                  
# Intit. variables
n_max_ecg = 3
n_max_eog = 2
title = 'Sources related to %s artifacts (red)'   

#%% Setup csv for logging removed components
dataCategories = ['subId','session','nRejct','nECG','nEOG']
saveFile = log_path+'/icaSummary.csv'          # Filename for save-data
csvWriter = csv.writer(open(saveFile, 'wt'), delimiter=';')            # The writer function to csv
csvWriter.writerow(dataCategories)

#%% Run
for ff in file_list: 
    print('Now reading file: '+ff)
    
    sub = ff[0:4]    
    condition = project_part
    session = ff[13]
    
    log = dict()
    log['subId'] = sub
    log['session'] = session

    subDir = path.join(meg_path,sub)
    ica_dir = op.join(subDir,'ica')
    out_fname = subDir+'/'+sub+'_RsEc'+session+'-ica_raw.fif'
    icaFname = project_part+'-ica.fif'
    
    chdir(subDir)

    if path.isfile(out_fname) and not overwrite_old_files:
        print('Will not overwrite existing file "'+out_fname+'". Pass!')
        continue

    in_fname = ff
    
#    if len(in_fname) > 1:
#        print('For '+ff[:10]+' there is '+str(len(in_fname))+' files to read!')
#        if path.isfile(out_fname):
#            print('File exist. Pass.')
#            continue
#        else:
#            print('Combining files to output: '+out_fname)
    
    #read ica from file
    if not path.isfile(op.join(ica_dir,icaFname)):
        print('ICA file '+icaFname+' does not exist. Pass!')
        continue
    else:
        print('found ICA file. Proceeding...')
    
    ica = read_ica(op.join(ica_dir,icaFname))
    ica.labels_ = dict()
    
    # LOAD DATA
#    for idx, ffs in enumerate(in_fname):
#        print(ffs)
#        if idx == 0:
    raw =  read_raw_fif(op.join(raw_path,ff), preload=True)
#        else:
#            raw.append(read_raw_fif(ffs, preload=True))

    picks_meg = pick_types(raw.info, meg=True, eeg=False, eog=False, emg=False, misc=False, 
                           stim=False, exclude='bads')
    picks_eXg = pick_types(raw.info, meg=False, eeg=False, eog=True, ecg = True, emg=False, misc=False, 
                           stim=False, exclude='bads')
    raw.filter(1, 40, n_jobs=3, picks=picks_eXg)
    raw.notch_filter(50, n_jobs=3, picks=picks_eXg)
    
    ecg_epochs = create_ecg_epochs(raw, ch_name='ECG003', tmin=-.5, tmax=.5)    #, picks=picks)
    ecg_inds, ecg_scores = ica.find_bads_ecg(ecg_epochs, method='ctps', verbose=False)
    
    ecg_scores_fig = ica.plot_scores(ecg_scores, exclude=ecg_inds, title=title % 'ecg', show=False)
    ecg_scores_fig.savefig(op.join(ica_dir,condition+session+'_ICA_ecg_comp_score.png'))

    if ecg_inds:
        show_picks = np.abs(ecg_scores).argsort()[::-1][:5]
        
        ecg_source_fig = ica.plot_sources(raw, show_picks, exclude=ecg_inds,
                         title=title % 'ECG', show=False)
        ecg_source_fig.savefig(op.join(ica_dir,condition+session+'_ICA_ecg_source.png'))

        ecg_comp_fig = ica.plot_components(ecg_inds, title=title % 'ecg', colorbar=True, show=False)
        ecg_comp_fig.savefig(op.join(ica_dir,condition+session+'_ICA_ecg_comp.png'))
        
    # estimate average artifact
    ecg_evoked = ecg_epochs.average()
    # plot ECG sources + selection
    ecg_evo_fig1 = ica.plot_sources(ecg_evoked, exclude=ecg_inds, show=False)    
    ecg_evo_fig1.savefig(op.join(ica_dir,condition+session+'_ICA_ecg_latSource.png'))
    
    ecg_inds = ecg_inds[:n_max_ecg]
    ica.exclude += ecg_inds
    log['nECG'] = len(ecg_inds)
    
    # Find EOG artifacts
    
    eog_inds, eog_scores = ica.find_bads_eog(raw)
    
    eog_scores_fig = ica.plot_scores(eog_scores, exclude=eog_inds, title=title % 'eog', show=False)
    eog_scores_fig.savefig(op.join(ica_dir,condition+session+'ICA_eog_comp_score.png'), show=False)

    if eog_inds:
        show_picks = np.abs(eog_scores[0]).argsort()[::-1][:5]
#        eog_source_fig = ica.plot_sources(raw, show_picks, exclude=eog_inds,
#                         title=title % 'EOG',show=False)
#                         
#        show_picks = np.abs(eog_scores[1]).argsort()[::-1][:5]
#        ica.plot_sources(raw, show_picks, exclude=eog_inds,
#                         title=title % 'EOG')
                         
        eog_comp_fig = ica.plot_components(eog_inds, title="Sources related to EOG artifacts",
                                colorbar=True,show=False)
        eog_comp_fig.savefig(op.join(ica_dir,condition+session+'ICA_eog_comp.png'))

    # estimate average artifact
    eog_evoked = create_eog_epochs(raw, tmin=-.75, tmax=.75, picks=picks_meg).average()
                                   
#    ica.plot_sources(eog_evoked, exclude=eog_inds)

    # plot EOG sources + selection
    eog_evo_fig = ica.plot_overlay(eog_evoked, exclude=eog_inds, show=False)  # plot EOG cleaning
    eog_evo_fig.savefig(op.join(ica_dir,condition+session+'ICA_eog_cleaning.png'))
        
    eog_inds = eog_inds[:n_max_eog]
    ica.exclude += eog_inds
    log['nEOG'] = len(eog_inds)
    
    # Apply the solution to Raw, Epochs or Evoked like this:
    log['nRejct'] = len(ica.exclude)
    
    
    raw_ica = ica.apply(raw)
    raw_ica.save(out_fname, overwrite=overwrite_old_files)
    
    csvWriter.writerow([log[category] for category in dataCategories])
    
    print('----------- DONE (Sub: '+sub+' Session: '+session+') -----------------')
    
print('----------- ALL DONE -----------------')

