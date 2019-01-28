#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Get beta-course from ROI
Created on Mon Jan 28 11:17:45 2019. @author: mikkel
"""
import numpy as np
from mne import find_events, read_labels_from_annot
import matplotlib.pyplot as plt
import os.path as op
from os import listdir
import sys
import mne
sys.path.append('/home/mikkel/PD_motor/tap/scripts/functions')
from sensorymotorROI import make_sensorymotorROI
from scipy.signal import hilbert
import scipy.io as sio


#%% Overwrite
overwrite = False

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
#subjects = ['0320','0313']

#%%
for sub in subjects:
    print('sub: '+sub)
    srcFile = op.join(src_dir,sub+'-ico4-src.fif')
    src = mne.read_source_spaces(srcFile)

    for con in ['RsEc1','RsEc2']:
        outmatfname = op.join(data_path,sub,sub+'_'+con+'-hilbt.mat')
        hilb = dict()
        
        for hemi in ['lh','rh']:
            fname = op.join(data_path,sub,sub+'_'+con+'-dSPM-lh.stc')
            stc = mne.read_source_estimate(fname)
            
            lab = make_sensorymotorROI(sub, subjects_dir, hemi=hemi) 
#
#            label_tcM = stc.extract_label_time_course(lab, src ,mode='mean')
#            label_tcM = np.float64(label_tcM)
#            label_tcMF = stc.extract_label_time_course(lab, src ,mode='mean_flip')[0,:]
            label_tcPCA = stc.extract_label_time_course(lab, src ,mode='pca_flip')[0,:]
            label_tcPCA  = np.float64(label_tcPCA)
#            label_tcX = stc.extract_label_time_course(lab, src ,mode='max')[0,:]
#            label_tcX  = np.float64(label_tcX)

#            t=np.arange(0,label_tcM.size)/1000.0
#            plt.plot(t,label_tcM)
#            plt.plot(t,label_tcMF)
#            plt.plot(t,label_tcPCA)
#            plt.plot(t,label_tcX)

#            lbftM = mne.filter.filter_data(label_tcM, 1000,10,30,method='iir')
            lbftPCA = mne.filter.filter_data(label_tcPCA, 1000, 10, 30, method='iir')
#            lbftX = mne.filter.filter_data(label_tcX, 1000,10,30,method='iir')


#            plt.plot(t,lbftM[0,:], 'r--')
#            plt.plot(t,lbftPCA, 'b--')
#            plt.plot(t,lbftX, 'k--')
            
#            analytic = hilbert(lbftM)
#            envelope = np.abs(analytic)6
            analyticPCA = hilbert(lbftPCA)
            envelopePCA = np.abs(analyticPCA)            
            
            hilb[hemi] = envelopePCA;
#            analyticX = hilbert(lbftX)
#            envelopeX = np.abs(analyticX)
            
#            plt.plot(t,envelope[0,:],'r-')
#            plt.plot(t,envelopePCA,'b-')
#            plt.plot(t,envelopeX,'k-')

#            plt.hist(envelope[0,:],50, alpha=0.5)
#            plt.hist(envelopePCA,50, alpha=0.5)
#            plt.hist(envelopeX,50, alpha=0.5)
            
        sio.savemat(outmatfname, dict(hilb_rh=hilb['rh'],hilb_lh=hilb['lh']))
        
#END