#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Make average source maps based on peaks in ROI timeseries.
Created on Tue Feb 26 14:05:23 2019. @author: mikkel
"""
import numpy as np
#import matplotlib.pyplot as plt
import os.path as op
from os import listdir
import sys
import mne
sys.path.append('/home/mikkel/PD_motor/tap/scripts/functions')
import scipy.io as sio
import pickle

#%% Overwrite
overwrite = True     # Overwrite  data if it already exists
do_plot   = False    # Plot the source maps?

#%% dirs
data_path       = '/home/mikkel/PD_motor/rest_ec/meg_data'         
subjects_dir    = '/home/mikkel/PD_motor/fs_subjects_dir'
src_dir         = '/home/mikkel/PD_motor/tap/mri'

subjects_meg = listdir(data_path)
subjects_mri = [f for f in listdir(subjects_dir) if f.startswith('0')]

subjects = list(set(subjects_meg).intersection(subjects_mri))
subjects.sort()

# Manual input for single subjects
#subjects = ['0320','0313']

#%% Run
Xall = dict()
Xall['RsEc1'] = np.zeros([len(subjects),5124])
Xall['RsEc2'] = np.zeros([len(subjects),5124])

for i, sub in enumerate(subjects):
    print('sub: '+sub)
    srcFile = op.join(src_dir,sub+'-ico4-src.fif')
    src = mne.read_source_spaces(srcFile)
    morph = mne.compute_source_morph(src, sub, 'fsaverage', subjects_dir=subjects_dir,spacing=4)

    for c, con in enumerate(['RsEc1','RsEc2']):       
        outfname = op.join(data_path,sub,sub+'_'+con+'-avg.stc')
        if op.exists(outfname) and not overwrite:
            continue
        
        # Load dSPM
        stcfname = op.join(data_path,sub,sub+'_'+con+'-dSPM-lh.stc')
        stc = mne.read_source_estimate(stcfname)
        
        # Load peak index
        pkfname = op.join(data_path,sub,'pkidx'+str(c+1)+'.mat')
        temp = sio.loadmat(pkfname)
        pkidx = (temp['maxidx'+str(c+1)])
        pkidx = pkidx-1                     # Correct for non-zero indexing
        
        # Select data and average
        tempdat = np.mean(stc.data[:,pkidx],1)
        mnstc = mne.SourceEstimate(tempdat, list([src[0]['vertno'],src[1]['vertno']]), tmin=0, tstep=1)
        
        if do_plot:
            mnstc.plot(sub,subjects_dir=subjects_dir,hemi='both', time_viewer=True)
        
        mnstc.save(outfname)
        
        # Morph
        mnstc_fsavg = morph.apply(mnstc, mri_resolution=False)

        Xall[con][i,:] = np.squeeze(mnstc_fsavg.data)
                  
# Save all data
with open('/home/mikkel/PD_motor/rest_ec/groupanalysis/Xall.pkl', 'wb') as f:
    pickle.dump(Xall, f, pickle.HIGHEST_PROTOCOL)  
print('DONE')

# END        