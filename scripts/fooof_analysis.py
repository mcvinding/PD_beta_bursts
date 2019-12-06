#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FOOOF analysis of powerspectra in PD-BB project.
Created on Wed Dec  4 16:19:05 2019. @author: mikkel
"""
# Import fooof functions for creating spectra and managing parameters
import numpy as np
import scipy.io
from fooof.sim.params import param_sampler, param_iter, Stepper
from fooof.sim.gen import gen_power_spectrum, gen_group_power_spectra
from fooof.analysis import get_band_peak, get_band_peak_group
import matplotlib.pyplot as plt

# Import some fooof plotting functions
from fooof.plts.spectra import plot_spectrum, plot_spectra
from fooof.plts.spectra import plot_spectrum_shading, plot_spectra_shading

from fooof import FOOOF
from fooof import FOOOFGroup

#%% SETTINGS
# Define frequency bands of interest
theta_band = [4, 8]
alpha_band = [8, 12]
beta_band = [13, 30]

groups = ['ctrl_dat1','ptns_dat1','ctrl_dat2','ptns_dat2']

#%% Read data
dat = scipy.io.loadmat('/home/mikkel/PD_motor/rest_ec/groupanalysis/PSD_data.mat')

#%% Do FOOOF
alpha_dict = dict()
beta_dict = dict()
theta_dict = dict()

for jj in groups:
    
    #Init.
    alpha_mat = np.zeros([len(dat[jj]),3])
    beta_mat = np.zeros([len(dat[jj]),3])
    theta_mat = np.zeros([len(dat[jj]),3])
    
    for ii in range(len(dat[jj])):
        print(ii)
        
        freq = dat[jj][ii][0][0][0][2][0,]
        psd = dat[jj][ii][0][0][0][3][0,]
    
        ## PLOT
        ##plot_spectrum(fs, psd[0,])
        #plot_spectrum(freq, psd_mat)
        #plot_spectrum_shading(freq, psd_mat, [13, 30], log_powers=True)
    
        # Initialize FOOOF object
        fm = FOOOF(max_n_peaks=6, min_peak_height=0.05, peak_width_limits = [.75, 20], aperiodic_mode='fixed')
    
        # Define frequency range across which to model the spectrum
        freq_range = [3, 35]
    
        # Fit model
        fm.fit(freq, psd, freq_range)
        
        # Model the power spectrum with FOOOF, and print out a report
        fm.print_results()
#        fm.plot()
    
        # Extract band-specific oscillations from the FOOOF model
        beta_mat[ii,] = get_band_peak(fm.peak_params_, beta_band)
        alpha_mat[ii,] = get_band_peak(fm.peak_params_, alpha_band)
        theta_mat[ii,] = get_band_peak(fm.peak_params_, theta_band)
    
    beta_dict[jj] = beta_mat.copy()
    alpha_dict[jj] = alpha_mat.copy()
    theta_dict[jj] = theta_mat.copy()


#%% Plot
for jj in groups:
    plt.figure(1)
    plt.hist(beta_dict[jj][:,1], 10); plt.title('Max peak pow')
    plt.figure(2)
    plt.hist(beta_dict[jj][:,0], 10); plt.title('Max peak freq')
    plt.figure(3)
    plt.hist(beta_dict[jj][:,2], 10); plt.title('Max peak width')

for jj in groups:
    plt.figure(4)
    plt.hist(theta_dict[jj][:,1], 10); plt.title('Max peak pow')
    plt.figure(5)
    plt.hist(theta_dict[jj][:,0], 10); plt.title('Max peak freq')
    plt.figure(6)
    plt.hist(theta_dict[jj][:,2], 10); plt.title('Max peak width')

for jj in groups:
    plt.figure(7)
    plt.hist(alpha_dict[jj][:,1], 10); plt.title('Max peak pow')
    plt.figure(8)
    plt.hist(alpha_dict[jj][:,0], 10); plt.title('Max peak freq')
    plt.figure(9)
    plt.hist(alpha_dict[jj][:,2], 10); plt.title('Max peak width')    
    
#%%     
# Create and save out a report summarizing the results across the group of power spectra
fg.save_report()

# Save out FOOOF results for further analysis later
fg.save(file_name='fooof_group_results', save_results=True)