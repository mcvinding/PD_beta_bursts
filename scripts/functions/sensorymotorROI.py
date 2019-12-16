#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Function for making sensory motor ROI.
Created on Fri Jan 18 14:24:35 2019. @author: mikkel
"""
import mne
from surfer import Brain

# Make sensory-motor ROI based on aparc.a2009s
def make_sensorymotorROI(subj, subjects_dir, hemi='both'):
    atlas = 'aparc.a2009s'    
    roilab = ['G_precentral', 'G_postcentral','S_central']

    labels = mne.read_labels_from_annot(subj, parc=atlas, subjects_dir=subjects_dir, hemi=hemi)
    ROI_labels = [lab for lab in labels if any(x in lab.name for x in roilab)]
    
    for i, lab in enumerate(ROI_labels):
        if i==0:
            ROI = lab
        else:
            ROI = ROI+lab
            
    ROI.name = 'sensorymotor-'+hemi    
    
    return ROI
    

## Plot ROI
#def plot_roi(subj, roi_label=None, subjects_dir, hemi='lh', fname=None):
#    if not roi_label:
#        roi_label = make_sensorymotorROI(subj, subjects_dir, hemi=hemi) 
#
#    brain = Brain('fsaverage', hemi, 'inflated', subjects_dir=subjects_dir,
#              cortex='bone', background='white', size=(800, 600))
#    brain.add_label(roi_label, alpha=1,borders=True)
#    brain.add_label(roi_label, alpha=0.5,borders=False)
#
#    if fname:
#        brain.save_single_image(fname)    
#

#END