#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Do MNE on resting state data. 
@author: mikkel
"""
import matplotlib
matplotlib.use('Agg')
from os import listdir
import os.path as op
from mne import read_forward_solution, read_bem_solution, read_source_spaces, read_trans, read_cov
from mne import make_bem_model, make_bem_solution, make_forward_solution
from mne import write_bem_solution, write_forward_solution
from mne.io import Raw
from mne.minimum_norm import make_inverse_operator, write_inverse_operator, apply_inverse_raw, read_inverse_operator
import time

#%% Overwrite
overwrite = True

#%% Paths, etc.
data_path       = '/home/mikkel/PD_motor/rest_ec/meg_data'         
subjects_dir    = '/home/mikkel/PD_motor/fs_subjects_dir'
src_dir         = '/home/mikkel/PD_motor/tap/mri'
trans_dir       = '/home/mikkel/PD_motor/tap/trans_files'  # I should be able to re-use as it is same didigizer data and same head cooridnate system
cov_dir         = '/home/mikkel/PD_motor/empty_room/noise_cov'

subjects_meg = listdir(data_path)
subjects_mri = [f for f in listdir(subjects_dir) if f.startswith('0')]

subjects = list(set(subjects_meg).intersection(subjects_mri))
subjects.sort()

# Manual input for single subjects
subjects = ['0406']

#%% Initiate values
startTrigger = 1
stopTrigger  = 64

conductivity = (0.3,)                               # for single layer
#lambda2 = 0.0000001

snr = 2.0                           # Standard assumption for average data is 3.0
lambda2 = 1.0 / snr ** 2

#%% SUBJECT LOOP HERE ####
for ii, sub in enumerate(subjects):
    print(sub)
    subDir = op.join(data_path,sub)
    rawfiles = [op.join(subDir,f) for f in listdir(subDir) if 'ica_raw.fif' in f]
    
    # Subject relevant files
    srcFile = op.join(src_dir,sub+'-ico4-src.fif')
    bemFile = op.join(subjects_dir,sub,'bem','bem-sol.fif')
    traFile = op.join(trans_dir,sub,sub+'-trans.fif' )                  #[!]
    covFile = op.join(cov_dir,sub,sub+'_combi-cov.fif')
    
    exisitngfiles = [op.join(subDir,f) for f in listdir(subDir) if 'dSPM' in f]
    if len(exisitngfiles) and not overwrite:
        print 'Output already exists. Will not overwrite!'
        continue
    
    ### Make bem (or read)
    if not op.isfile(bemFile) or overwrite:
        model = make_bem_model(subject=sub, ico=4, conductivity=conductivity, subjects_dir=subjects_dir)
        bem = make_bem_solution(model)
        write_bem_solution(op.join(subjects_dir,sub,'bem','bem-sol.fif'),bem)
        print('Wrote BEM file '+str(op.join(subjects_dir,sub,'bem','bem-sol.fif')))
    else:
        bem = read_bem_solution(bemFile)

    ### Read source space
    src = read_source_spaces(srcFile)
    
    ### read trans
    trans = read_trans(traFile)
    
    ### Read noise covariance from empty room recordings
    noise_cov = read_cov(covFile)    
    
    #### LOOP THROUGH RAWFILES
    for ff, rawFile in enumerate(rawfiles):
        print(rawFile)
        
        # Initiate filenames for loading and saving    
        fwdFile = rawFile[:-12]+'-fwd.fif'
        invFile = rawFile[:-12]+'-inv.fif'
    
        outFile_raw  = rawFile[:-11]+'crop-raw.fif'
        outFile_dSPM = rawFile[:-11]+'dSPM'
        outFile_MNE  = rawFile[:-11]+'MNE'                             

        # Read raw data
        if not op.isfile(outFile_raw) or overwrite:
            raw = Raw(rawFile, preload=True)
            eve = find_events(raw)
            
            try:
                startSam = eve[eve[:,2] == startTrigger,0][0]
            except:
                startSam = raw.first_samp

            try:
                stopSam = eve[eve[:,2] == stopTrigger,0][0]
            except:
                stopSam = raw.last_samp 

#            if not startSam:
#                startSam = raw.first_samp
#            if not stopSam:
#                stopSam = raw.last_samp 
            
            raw_crop = raw.crop(tmin=(startSam - raw.first_samp ) / raw.info['sfreq'],
                                tmax=(stopSam - raw.first_samp ) / raw.info['sfreq'])
            
#            raw_crop.filter(14,25)
            
            raw_crop.save(outFile_raw, overwrite=True)
            
            del raw             # Clear memory
            
        raw_crop =  Raw(outFile_raw, preload=False)
    
        # Make forward model
        if not op.isfile(fwdFile) or overwrite:
            fwd = make_forward_solution(raw_crop.info, trans=trans, src=src, bem=bem,
                                        meg=True, eeg=False,
                                        mindist=0.0, n_jobs=3)
            write_forward_solution(fwdFile, fwd, overwrite=True)
        else:
            fwd = read_forward_solution(fwdFile)

        # Make (or read) inverse operator
        if not op.isfile(invFile) or overwrite:
            inv = make_inverse_operator(raw_crop.info, fwd, noise_cov)
            write_inverse_operator(invFile,inv)
        else:
            inv = read_inverse_operator(invFile)
            
        # Do source recon
        t0 = time.time()
        stc_dSPM    = apply_inverse_raw(raw_crop, inv, lambda2, method="dSPM")
#        stc_MNE     = apply_inverse_raw(raw_crop, inv, lambda2, method="MNE")
        dt = time.time() - t0
        print('Time elapsed: '+str(dt/60.0)+' min')
        

#       Save
        stc_dSPM.save(outFile_dSPM)
#        stc_MNE.save(outFile_MNE)

