#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FOOOF analysis of powerspectra in PD-BB project.
Created on Wed Dec  4 16:19:05 2019. @author: mikkel
"""
# Import fooof functions for creating spectra and managing parameters
import numpy as np
import scipy.io
import pandas
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
aper_dict = dict()
rsq_dict = dict()

for jj in groups:
    
    #Init.
    alpha_mat = np.zeros([len(dat[jj]),3])
    beta_mat = np.zeros([len(dat[jj]),3])
    theta_mat = np.zeros([len(dat[jj]),3])
    aper_mat = np.zeros([len(dat[jj]),2])
    
    rsq_mat = np.zeros([len(dat[jj]),1])
    
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
        freq_range = [1, 45]
    
        # Fit model
        fm.fit(freq, psd, freq_range)
        
        # Model the power spectrum with FOOOF, and print out a report
#        fm.print_results()
#        fm.plot()
    
        # Extract band-specific oscillations from the FOOOF model
        beta_mat[ii,] = get_band_peak(fm.peak_params_, beta_band)
        alpha_mat[ii,] = get_band_peak(fm.peak_params_, alpha_band)
        theta_mat[ii,] = get_band_peak(fm.peak_params_, theta_band)
        aper_mat[ii,] = fm.get_params('aperiodic_params')
        rsq_mat[ii] = fm.get_params('r_squared')
        
    beta_dict[jj] = beta_mat.copy()
    alpha_dict[jj] = alpha_mat.copy()
    theta_dict[jj] = theta_mat.copy()
    aper_dict[jj] = aper_mat.copy()
    rsq_dict[jj] = rsq_mat.copy()

#%% Arrange dataframe and export

df_beta = pandas.DataFrame()
df_alpha = pandas.DataFrame()
df_theta = pandas.DataFrame()

for con in beta_dict:
    temp = beta_dict[con].copy()
    print(con)
    
    session = np.ones(len(temp)) if con[-1]=='1' else np.ones(len(temp))*2
    group = np.ones(len(temp)) if con[:4]=='ptns' else np.ones(len(temp))*2
    
    d_temp = {'peak_freq': temp[:,0],
              'peak_pow': temp[:,1],
              'peak_width': temp[:,2],
              'group': group,
              'session': session}
    db_temp = pandas.DataFrame(d_temp)
    df_beta = df_beta.append(db_temp)

for con in alpha_dict:
    temp = alpha_dict[con].copy()
    print(con)
    
    session = np.ones(len(temp)) if con[-1]=='1' else np.ones(len(temp))*2
    group = np.ones(len(temp)) if con[:4]=='ptns' else np.ones(len(temp))*2
    
    d_temp = {'peak_freq': temp[:,0],
              'peak_pow': temp[:,1],
              'peak_width': temp[:,2],
              'group': group,
              'session': session}
    db_temp = pandas.DataFrame(d_temp)
    df_alpha = df_alpha.append(db_temp)
    
    
for con in theta_dict:
    temp = theta_dict[con].copy()
    print(con)
    
    session = np.ones(len(temp)) if con[-1]=='1' else np.ones(len(temp))*2
    group = np.ones(len(temp)) if con[:4]=='ptns' else np.ones(len(temp))*2
    
    d_temp = {'peak_freq': temp[:,0],
              'peak_pow': temp[:,1],
              'peak_width': temp[:,2],
              'group': group,
              'session': session}
    db_temp = pandas.DataFrame(d_temp)
    df_theta = df_theta.append(db_temp)

# Export    
df_beta.to_csv('/home/mikkel/PD_motor/rest_ec/groupanalysis/df_beta.csv', index=False, sep=';')
df_alpha.to_csv('/home/mikkel/PD_motor/rest_ec/groupanalysis/df_alpha.csv', index=False, sep=';')
df_theta.to_csv('/home/mikkel/PD_motor/rest_ec/groupanalysis/df_theta.csv', index=False, sep=';')

#%% Inspect
for p in rsq_dict:
    plt.figure(1)
    plt.hist(rsq_dict[p], alpha=.5, label=p); plt.title('R-squared')
plt.legend()    

for p in aper_dict:
    plt.figure(2)
    plt.hist(aper_dict[p][:,0], alpha=.5, edgecolor='k', label=p); plt.title('Intercept')
plt.legend()    

for p in aper_dict:    
    plt.figure(3)
    plt.hist(aper_dict[p][:,1], alpha=.5, edgecolor='k', label=p); plt.title('slope')
plt.legend()    

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Below here still in progress: move to somewhere else 

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
    

#%% Combined figure
colors = [[0,.5,1], [0,0,1], [1,.5,0], [1,0,0] ]
dat_b = [beta_dict['ctrl_dat1'][:,0], beta_dict['ctrl_dat2'][:,0], beta_dict['ptns_dat1'][:,0], beta_dict['ptns_dat2'][:,0]]
dat_a = [alpha_dict['ctrl_dat1'][:,0], alpha_dict['ctrl_dat2'][:,0], alpha_dict['ptns_dat1'][:,0], alpha_dict['ptns_dat2'][:,0]]
dat_t = [theta_dict['ctrl_dat1'][:,0], theta_dict['ctrl_dat2'][:,0], theta_dict['ptns_dat1'][:,0], theta_dict['ptns_dat2'][:,0]]

bins = np.linspace(4, 31, 28)    

plt.figure()
plt.hist(dat_b, bins, label='x', edgecolor='k', color=colors, alpha=1, ls='-')
plt.hist(dat_a, bins, label='x', edgecolor='k', color=colors, alpha=1, ls='--')
plt.hist(dat_t, bins, label='x', edgecolor='k', color=colors, alpha=1, ls=':')

#%%
plt.hist(beta_dict['ptns_dat1'][:,0], bins, alpha=0.5, label='x', edgecolor='k',hatch='/', facecolor='b', histtype='stacked')
plt.hist(beta_dict['ptns_dat2'][:,0], bins, alpha=0.5, label='x', edgecolor='k', facecolor='b', histtype='stacked')
plt.hist(beta_dict['ctrl_dat1'][:,0], bins, alpha=0.5, label='x', edgecolor='k',hatch='/', facecolor='r', histtype='stacked')
plt.hist(beta_dict['ctrl_dat2'][:,0], bins, alpha=0.5, label='x', edgecolor='k', facecolor='r', histtype='stacked')



plt.hist(theta_dict['ctrl_dat1'][:,0], bins, alpha=0.5, label='x', edgecolor='b')
plt.hist(alpha_dict['ctrl_dat1'][:,0], bins, alpha=0.5, label='x', edgecolor='b')
plt.hist(beta_dict['ctrl_dat2'][:,0], bins, alpha=0.5, label='x', edgecolor='b')
plt.hist(theta_dict['ctrl_dat2'][:,0], bins, alpha=0.5, label='x', edgecolor='b')
plt.hist(alpha_dict['ctrl_dat2'][:,0], bins, alpha=0.5, label='x', edgecolor='b')

pyplot.hist(y, bins, alpha=0.5, label='y')
pyplot.legend(loc='upper right')
pyplot.show()

    
    
    
#%%     
# Create and save out a report summarizing the results across the group of power spectra
fg.save_report()

# Save out FOOOF results for further analysis later
fg.save(file_name='fooof_group_results', save_results=True)