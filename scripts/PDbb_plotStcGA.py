#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Plot group level source maps.
Created on Tue Feb 26 16:43:42 2019. @author: mikkel
"""
import pickle
import mne
from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import mayavi

#%% Load data
with open('/home/mikkel/PD_motor/rest_ec/groupanalysis/Xall.pkl', 'wb') as f:
    Xall = pickle.load(f)  

#%% Paths
data_path       = '/home/mikkel/PD_motor/rest_ec/meg_data'         
subjects_dir    = '/home/mikkel/PD_motor/fs_subjects_dir'
src_dir         = '/home/mikkel/PD_motor/tap/mri'

subjects_meg = listdir(data_path)
subjects_mri = [f for f in listdir(subjects_dir) if f.startswith('0')]

subjects = list(set(subjects_meg).intersection(subjects_mri))
subjects.sort()

# Manual input for single subjects
#subjects = ['0320','0313']

# Seperate subject groups
map_controls=['0320','0328','0343','0360','0361','0362','0363','0364','0367','0371','0373','0374','0375','0376','0386','0389','0393','0398','0399','0403']
map_patients=['0313','0314','0319','0322','0325','0327','0332','0333','0339','0340','0342','0345','0352','0353','0355','0366','0377','0392','0397','0406']

ctrl_subs = set([s for s in map_controls if s in subjects])
ptns_subs = set([s for s in map_patients if s in subjects])

ctrl_idx =  [i for i, val in enumerate(subjects) if val in ctrl_subs]
ptns_idx =  [i for i, val in enumerate(subjects) if val in ptns_subs]

#%% Split data
Xptns1 = Xall['RsEc1'][ptns_idx,:]
Xctrl1 = Xall['RsEc1'][ctrl_idx,:]
Xptns2 = Xall['RsEc2'][ptns_idx,:]
Xctrl2 = Xall['RsEc2'][ctrl_idx,:]

# Average
src = mne.read_source_spaces('/home/mikkel/PD_motor/tap/mri/fsaverage-ico4-src.fif')
GAptns1 = mne.SourceEstimate(np.mean(Xptns1,0), list([src[0]['vertno'],src[1]['vertno']]), tmin=0, tstep=1)
GActrl1 = mne.SourceEstimate(np.mean(Xctrl1,0), list([src[0]['vertno'],src[1]['vertno']]), tmin=0, tstep=1)
GAptns2 = mne.SourceEstimate(np.mean(Xptns2,0), list([src[0]['vertno'],src[1]['vertno']]), tmin=0, tstep=1)
GActrl2 = mne.SourceEstimate(np.mean(Xctrl2,0), list([src[0]['vertno'],src[1]['vertno']]), tmin=0, tstep=1)

fig1 = mayavi.mlab.figure()
fig2 = mayavi.mlab.figure()
fig3 = mayavi.mlab.figure()
fig4 = mayavi.mlab.figure()

GAptns1.plot('fsaverage',subjects_dir=subjects_dir, surface='inflated', hemi='split', figure=list([fig1,fig2]))
GAptns2.plot('fsaverage',subjects_dir=subjects_dir, surface='inflated', hemi='split', figure=list([fig3,fig4]))
