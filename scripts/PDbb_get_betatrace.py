#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Get beta-course from ROI
Created on Mon Jan 28 11:17:45 2019. @author: mikkel
"""
import numpy as np
import os.path as op
from os import listdir
import sys
import mne
sys.path.append('/home/mikkel/PD_motor/rest_ec/scripts/functions')
from sensorymotorROI import make_sensorymotorROI
from sensorymotorROI2 import make_sensorymotorROI2

from scipy.signal import hilbert
import scipy.io as sio

#%% Overwrite
overwrite = False

#%% Paths, etc.
data_path       = '/home/mikkel/PD_motor/rest_ec/meg_data'         
subjects_dir    = '/home/mikkel/PD_motor/fs_subjects_dir'
src_dir         = '/home/mikkel/PD_motor/tap/mri'
trans_dir       = '/home/mikkel/PD_motor/tap/trans_files'  # I should be able to re-use as it is same digitzer data and same head cooridnate system
cov_dir         = '/home/mikkel/PD_motor/empty_room/noise_cov'

subjects_meg = listdir(data_path)
subjects_mri = [f for f in listdir(subjects_dir) if f.startswith('0')]

subjects = list(set(subjects_meg).intersection(subjects_mri))
subjects.sort()

# Manual input for single subjects
#subjects = ['0320','0313']

#%% Run
for sub in subjects:
    print('sub: '+sub)
    srcFile = op.join(src_dir,sub+'-ico4-src.fif')
    src = mne.read_source_spaces(srcFile)

    for con in ['RsEc1','RsEc2']:
        outhilbt = op.join(data_path,sub,sub+'_'+con+'-hilbt.mat')      # Hilbert envelope
        outrawtc = op.join(data_path,sub,sub+'_'+con+'-rawtc.mat')      # Raw time-series
        outrawft = op.join(data_path,sub,sub+'_'+con+'-rawft.mat')      # Band-pass filtered time-series

        hilb = dict()
        rawtc = dict()
        filtc = dict()
        
        fname = op.join(data_path,sub,sub+'_'+con+'-dSPM-lh.stc')
        stc = mne.read_source_estimate(fname)
        
        for hemi in ['lh','rh']:

            lab = make_sensorymotorROI(sub, subjects_dir, hemi=hemi)

            # Save label
            lab.save(op.join(subjects_dir, sub, 'label',hemi+'.sensmotor.label'))

            # Extract label time-series
            label_tc = stc.extract_label_time_course(lab, src ,mode='pca_flip')[0,:]
            label_tc  = np.float64(label_tc)
            
            lbft = mne.filter.filter_data(label_tc, 1000, 13, 30, method='fir',n_jobs=3)

            analytic = hilbert(lbft)
            envelope = np.abs(analytic)         
            
            hilb[hemi] = envelope
            rawtc[hemi] = label_tc
            filtc[hemi] = lbft
            
        if not op.exists(outhilbt) or overwrite:    
            sio.savemat(outhilbt, dict(hilb_rh=hilb['rh'],hilb_lh=hilb['lh']))
        if not op.exists(outrawtc) or overwrite:
            sio.savemat(outrawtc, dict(raw_rh=rawtc['rh'],raw_lh=rawtc['lh']))
        if not op.exists(outrawft) or overwrite:
            sio.savemat(outrawft, dict(rwft_rh=filtc['rh'],rwft_lh=filtc['lh']))
        
#END