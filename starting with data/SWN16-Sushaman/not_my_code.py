#!/usr/bin/env python
# coding: utf-8

import sys, os
os.environ['PYSYN_CDBS']='/Users/iair/data/cdbs/'

import astropy
from astropy.table import Table
import numpy as np
import pandas as pd
import pysynphot as S
import random

from scipy import constants
import time as time1
import datetime

import matplotlib.pyplot as plt
#import matplotlib.lines as mlines
#import matplotlib.patches as mpatches

import emcee
# import beautifulcorner

# Useful Constants

k_B = constants.value('Boltzmann constant')  # in [J]*[K]^(-1)
h = constants.value('Planck constant') # in [kg]*[m]^2*[s]^(-1)
c = astropy.constants.c.value # in [m]*[s]^(-1)
m = 15

distance = 26.4e6 #in pc units
distance_cm = 8.14618882e25 # in cm
redshift = 0.00494
t0 = 2457651.69028 # VSNET Discovery JD (actual first image).

# Paths

swift_bp_path = '/Users/iair/data/filter_curves/UVOTfilters_Breeveld2011/'
atlas_bp_path = '/Users/iair/data/filter_curves/ATLAS/'
data_path = '/Users/iair/Documents/TAU/Research/SN 2016gkg/OrMoalem/'
mcmc_path = '/Users/iair/Documents/TAU/Research/SN 2016gkg/OrMoalem/mcmc_results/'

# Load the Data

df_I = pd.read_csv(data_path + 'fig1data_extinction_corrected.csv')

#separate data with or without lower bound
df_I_limited = df_I[df_I['Is_Limit'] > 0.5] #limited values get Is_Limit value 1
df_I_unlimited = df_I[df_I['Is_Limit'] < 0.5] #unlimited values get Is_Limit value 0
df_I_unlimited_copy = df_I_unlimited #for future usage 

#set 'Filter' column to be the index of the dataframes
df_I_limited = df_I_limited.set_index("Filter")
df_I_unlimited = df_I_unlimited.set_index("Filter")

#load amateur data
df_A = pd.read_csv(data_path + 'input_amateur_start_also_averaged.csv')
#print(df_A)

#separate data with or without lower bound
df_A_limited = df_A[df_A['Is_Limit'] > 0.5] #limited values get Is_Limit value 1
df_A_unlimited = df_A[df_A['Is_Limit'] < 0.5] #unlimited values get Is_Limit value 0

# Order the Data 

filters = ['V', 'o', 'UVW1', 'U', 'B', 'UVW2', 'UVM2', 'g', 'r', 'i']
colors = ['b', 'c', 'm', 'y', 'k', 'g', 'slateblue', 'violet', 'olive', 'lightgreen']

bandpasses = {}
for filtname in filters:
    if filtname in ['g','r','i']: 
        newfiltname = 'sdss,'+filtname
    else:
        newfiltname = filtname
    if newfiltname in ['UVW1','UVW2','UVM2']:
        bandpasses[filtname] = S.FileBandpass(swift_bp_path + newfiltname+'_UVOT.txt')
    elif newfiltname == 'o':
        bandpasses[filtname] = S.FileBandpass(atlas_bp_path + 'orange_in_aa.dat')
    elif newfiltname == 'Clear':
        bandpasses[filtname] = S.ObsBandpass('sdss,r')
    else:
        bandpasses[filtname] = S.ObsBandpass(newfiltname)

df_time_mcmc = pd.DataFrame()
df_mag_mcmc = pd.DataFrame()
df_magErr_mcmc = pd.DataFrame()

for f in filters :
    df_time_mcmc.loc[:,f] = ((df_I_unlimited.JD.loc[f]-df_I.JD.loc[0])*24*3600/(1+redshift)).reset_index(drop=True) #in seconds
    df_mag_mcmc.loc[:,f] = (df_I_unlimited.Magnitude.loc[f]).reset_index(drop=True)
    df_magErr_mcmc.loc[:,f] = (df_I_unlimited.Magnitude_Error.loc[f]).reset_index(drop=True)

#amateur's data is in filter V
df_time_mcmc.loc[:,'V'].append((df_A_unlimited.JD-df_I.JD.loc[0])*24*3600/(1+redshift)).reset_index(drop=True)
df_mag_mcmc.loc[:,'V'].append(df_A_unlimited.Magnitude).reset_index(drop=True)
df_magErr_mcmc.loc[:,'V'].append(df_A_unlimited.Magnitude_Error).reset_index(drop=True)

data_t = {}
data_mag = {}
data_magerr = {}

for f in filters:
    data_t[f] = df_time_mcmc[f].dropna().to_numpy()
    data_mag[f] = df_mag_mcmc[f].dropna().to_numpy()
    data_magerr[f] = df_magErr_mcmc[f].dropna().to_numpy()

# Variables for mcmc simulation
R_star = 500 # radius of supernove in unit of R_sun
M_ej = 15 # ejected mass in units if M_sun
E_exp = 1 # explosion energy in unit of 10^51 [erg]
toffset_d = 0.01 # offset in days relative to discovery (t0)

M15 = M_ej/15 
R500 = R_star/500 
E51 = E_exp 

# special times for model
t0 = 155*(M15**0.23)*(R500**1.39)*(E51**(-0.81)) - toffset_d*24*3600 # formulas according to the article + offset
ts = 3.6*3600*(M15**0.44)*(R500**1.49)*(E51**(-0.56)) - toffset_d*24*3600
tc = 3.2*24*3600*(M15**0.97)*(R500**2.02)*(E51**(-1.19)) - toffset_d*24*3600
trec = 16.6*24*3600*(M15**0.22)*(R500**0.76)*(E51**(-0.43)) - toffset_d*24*3600
t_Rc = (R_star*astropy.constants.R_sun.value)/astropy.constants.c.value # in units of seconds

max_days_to_fit = 3

# Additional variables for future usage
step_lambda = 100
lambda_min = 1000 
lambda_max = 12000

n = 50

##############################################################################

filtclr={}
filtclr['U']='#3c0072'
filtclr['B']='#0057ff'
filtclr['V']='#6DE600'
filtclr['R']='#ff7000'
filtclr['I']='#80000d'
filtclr['Bessell-B']='#0057ff'
filtclr['Bessell-V']='#6DE600'
filtclr['Bessell-R']='#ff7000'
filtclr['Bessell-I']='#80000d'
filtclr['up']='#2c0080'
filtclr['gp']='#00ccff'
filtclr['rp']='#ff7d00'
filtclr['ip']='#90002c'
filtclr['zs']='black'
filtclr['u']='#2c0080'
filtclr['g']='#00ccff'
filtclr['r']='#ff7d00'
filtclr['i']='#90002c'
filtclr['z']='black'
filtclr['w']='#a883b8'
filtclr['SDSS-U']='#2c0080'
filtclr['SDSS-G']='#00ccff'
filtclr['SDSS-R']='#ff7d00'
filtclr['SDSS-I']='#90002c' 
filtclr['Pan-Starrs-Z']='black'
filtclr['Pan-Starrs']='black'
filtclr['H-Alpha']='grey'
filtclr['air']='black'
filtclr['Clear']='white'
filtclr['G']='#a883b8'
filtclr['o']='y'

filtclr['UVW2'] = 'pink'
filtclr['UVM2'] = 'hotpink'
filtclr['UVW1'] = 'deeppink'

# Modified:
filtclr['V']='#46CC00'
filtclr['R']='#CC5900'
filtclr['z']='grey' 
filtclr['I']='grey'
    
##############################################################################

def cal_L2(t, M15, R500, E51, t0, ts, trec):
    if t <= 0 :
        L = 1.8*(10**45)*(M15**(-0.65))*(R500**(-0.11))*(E51**1.37)*np.exp(-0.35*((t/t0)**2)+0.15*(t/t0))
    if (t > 0 and t <= t0) :
        L = 1.8*(10**45)*(M15**(-0.65))*(R500**(-0.11))*(E51**1.37) #t<<t0
    if (t > t0 and t <= ts) : 
        t_hr = t/3600
        L = 2.7*(10**43)*(M15**(-0.34))*(R500**(1.74))*(E51**0.29)*(t_hr**(-4/3)) #t0<<t<<ts
    if (t > ts and t <= trec) :
        t_dy = t/(3600*24)
        L = 1.6*(10**42)*(M15**(-0.78))*(R500**(0.28))*(E51**0.84)*(t_dy**(-0.35)) #ts<<t<<trec
    if t > trec :
        return
    return L

##############################################################################

def cal_temp2(t, M15, R500, E51, t0, ts, tc, trec):
    if t <= t0 :
        T = 4.3*(10**5)*(M15**(-0.17))*(R500**(-0.52))*(E51**0.35) #t<<t0
    if (t > t0 and t <= ts) : 
        t_hr = t/3600
        T = (10**5)*(M15**(-0.07))*(R500**(0.1))*(E51**(-0.01))*(t_hr**(-0.45)) #t0<<t<<ts
    if (t > ts and t <= tc) :
        t_dy = t/(3600*24)
        T = 3*(10**4)*(M15**(-0.11))*(R500**(-0.04))*(E51**(0.04))*(t_dy**(-0.35)) #ts<<t<<tc
    if (t > tc and t <= trec) :
        t_dy = t/(3600*24)
        T = 4.1*(10**4)*(M15**(0.13))*(R500**(0.46))*(E51**(-0.25))*(t_dy**(-0.6)) #tc<<t<<trec
    if t > trec :
        return
    return T

##############################################################################

def calc_mag(step_lambda, lambda_min, lambda_max, temperature, lbol, f):
    
    #Create dataframe with frequency and flux column 
    #df = pd.DataFrame(np.arange(lambda_min, lambda_max, step_lambda), columns=['lambda'])
    df = {'lambda': np.arange(lambda_min, lambda_max, step_lambda)}
    df['nu'] = ((10**10)*c)/df['lambda'] #in 1/s
    df['F_lambda'] = np.tile(np.nan,len(df['lambda'])) # create an empty column
    
    #print("5a: %s seconds ---" % (time1.time() - start_time))
    for i in range(0, len(df['lambda'])): #loop for frequency
        temperature_col = temperature*((1+((h*df['nu'][i])/((3*k_B*temperature)**(-0.2*m))))**(-1/m))
        L_nu = 0.9*lbol*(15/(np.pi**4))*((h/(k_B*temperature_col))**4)*(df['nu'][i]**3)*((np.exp((h*df['nu'][i])/(k_B*temperature_col))-1)**(-1))
        L_lambda = (df['nu'][i]**2/c)*L_nu
        df['F_lambda'][i] = L_lambda/(4*np.pi*(distance_cm**2)*(10**10)) # flux in units of erg/(s*cm^2*Angstrom)
    #print("5b: %s seconds ---" % (time1.time() - start_time))
    
    sp = S.ArraySpectrum(df['lambda'], df['F_lambda'], fluxunits='flam') # flux in units of erg/(s*cm^2*Angstrom), wavelength in Angstrum  
    bp = bandpasses[f]
    if f in ['U','B','V','UVM2','UVW1','UVW2']:
        magsystem = 'VegaMag'
    elif f in ['g','r','i','o']:
        magsystem = 'ABMag'
    else:
        magsystem = 'ABMag'
    obs = S.Observation(sp,bp)
    try:
        mag = obs.effstim(magsystem)  
        #print("5c: %s seconds ---" % (time1.time() - start_time))
        return mag
    except: # For some parameters the flux is <0?
        return np.inf

##############################################################################

def correct_light_travel_time(f, t, t_Rc, R_star, n, step_lambda, lambda_min, lambda_max, expected_mag_before_correction):
    #only correct if t < 5*t_Rc
    if t < 5*t_Rc :   
        t_array = np.arange(t - t_Rc, t, t_Rc/n)
        integrand = []
        for t_i in t_array :
            lbol = cal_L2(t_i, M15, R500, E51, t0, ts, trec)
            temperature = cal_temp2(t_i, M15, R500, E51, t0, ts, tc, trec)
            mag = calc_mag(step_lambda, lambda_min, lambda_max, temperature, lbol, f)
            if mag == np.inf:
                return np.inf
            flux = 10**(-0.4*mag)
            integrand.append(flux*(1-(c*(t-t_i))/(R_star*astropy.constants.R_sun.value)))
        flux_corrected = (1/t_Rc)*np.trapz(integrand, t_array)
        mag_corrected = -2.5*np.log10(flux_corrected)
        return mag_corrected
    else:
        return expected_mag_before_correction

##############################################################################

def log_likelihood(theta, filters, data_t , data_mag, data_magerr):

    M15, R500, E51, toffset_d = theta
        
    if M15 < 0.1 or M15 > 2 or R500 < 0.1 or R500 > 4 or E51 < 0.1 or E51 > 5 or toffset_d < -1 or toffset_d > 1 :
        return -np.inf    
       
    chi2 = []
    normalization = []
    expected_mags = []
    
    #print("1: %s seconds ---" % (time1.time() - start_time))
    #print(filters)
    # Calculate expected magnitude
    for f in filters:
        t = data_t[f]
        mag = data_mag[f]
        magerr = data_magerr[f]
        #print(f)
        for i in range(0, len(t)): #Loop for time
            #print("2: %s seconds ---" % (time1.time() - start_time))
            lbol = cal_L2(t[i]-toffset_d*24*3600, M15, R500, E51, t0, ts, trec)
            #print("3: %s seconds ---" % (time1.time() - start_time))
            temperature = cal_temp2(t[i]-toffset_d*24*3600, M15, R500, E51, t0, ts, tc, trec)
            #print("4: %s seconds ---" % (time1.time() - start_time))
            expected_mag_before_correction = calc_mag(step_lambda, lambda_min, lambda_max, temperature, lbol, f)    
            if expected_mag_before_correction == np.inf:
                return -np.inf
            #print("5: %s seconds ---" % (time1.time() - start_time))
            
            # Correction taking in account the light travel time
            expected_mag = correct_light_travel_time(f, t[i]-toffset_d*24*3600, t_Rc, R_star, n, step_lambda, lambda_min, lambda_max, expected_mag_before_correction)
            if expected_mag == np.inf:
                return -np.inf
            #print("6: %s seconds ---" % (time1.time() - start_time))
            
            # Convert measured and expected magnitude to measured and expected flux for chi2 measurement:
            measured_flux = 10**(-0.4*mag[i]) 
            measured_fluxerr = abs(measured_flux*(-0.4)*np.log(10)*magerr[i]) 
            expected_flux = 10**(-0.4*expected_mag)

            chi2.append(((expected_flux-measured_flux)/measured_fluxerr)**2)
            normalization.append(np.log(2 * np.pi * measured_fluxerr**2))
            
            #print(" ------------- ")
            
    #time1.sleep(2)
    #print(" ------------- ")
    #print(" ------------- ")
    return -0.5 * np.sum(normalization + chi2)

##############################################################################

def log_flat_prior(theta):
    M15, R500, E51, toffset_d = theta

    if M15 > 0.6 or M15 < 2 or R500 > 0.1 or R500 < 4 or E51 > 0.1 or E51 < 5 or toffset_d > -1 or toffset_d < 1 :
        return 0.0
    else:
        return -np.inf 
    
def log_logflat_prior(theta):
    M15, R500, E51, toffset_d = theta

    if M15 > 0.6 or M15 < 2 or R500 > 0.1 or R500 < 4 or E51 > 0.1 or E51 < 5 or toffset_d > -1 or toffset_d < 1 :
        return (-np.log10(M15)) + (-np.log10(R500)) + (-np.log10(E51))
    else:
        return -np.inf 

##############################################################################

def log_posterior(theta, filters, df_time_mcmc , df_mag_mcmc, df_magErr_mcmc):
    return log_flat_prior(theta) + log_likelihood(theta, filters, df_time_mcmc , df_mag_mcmc, df_magErr_mcmc)

##############################################################################

def calc_model_mags(M15, R500, E51, toffset_d, filters, t):
        
    if M15 < 0.1 or M15 > 2 or R500 < 0.1 or R500 > 4 or E51 < 0.1 or E51 > 5 or toffset_d < -1 or toffset_d > 1 :
        return None
       
    expected_mags = {}
    
    for f in filters:
        expected_mags[f] = []
        for i in range(0, len(t)): #Loop for time
            lbol = cal_L2(t[i]-toffset_d*24*3600, M15, R500, E51, t0, ts, trec)
            temperature = cal_temp2(t[i]-toffset_d*24*3600, M15, R500, E51, t0, ts, tc, trec)
            expected_mag_before_correction = calc_mag(step_lambda, lambda_min, lambda_max, temperature, lbol, f)    
            if expected_mag_before_correction == np.inf:
                return None
            expected_mag = correct_light_travel_time(f, t[i]-toffset_d*24*3600, t_Rc, R_star, n, step_lambda, lambda_min, lambda_max, expected_mag_before_correction)
            if expected_mag == np.inf:
                return None
            expected_mags[f].append(expected_mag)
    return expected_mags
    
##############################################################################

nw = datetime.datetime.now()
nw_str = str(nw).replace('-','').replace(' ','_').replace(':','').split('.')[0]
summary_filename = mcmc_path + 'results_' + nw_str + '.txt'
chain_filename = mcmc_path + 'chain_' + nw_str + '.txt'

domcmc = True
loadfromfile = False
filetoload = 'mcmc_results/chain_20200219_140559.txt'

##############################################################################

if domcmc:
    
    data_t_for_fit = {}
    data_mag_for_fit = {}
    data_magerr_for_fit = {}
    for f in filters:
        w = data_t[f] <= max_days_to_fit*24*3600
        data_t_for_fit[f] = data_t[f][w]
        data_mag_for_fit[f] = data_mag[f][w]
        data_magerr_for_fit[f] = data_magerr[f][w]
      
    ndim = 4
    nwalkers = 8 
        
    # initialize walkers [M15, R500, E51, toffset_d]
    
    initial_guess = [1, 2, 2.5, 0.01]
    initial_spread = [0.2, 0.4, 0.5, 0.2] #לשחק אם הם לא מתכנסים
    starting_positions = np.random.randn(nwalkers, ndim)*initial_spread+initial_guess
    for wlkr in range(nwalkers):
        starting_positions[wlkr][3] = np.min([starting_positions[wlkr][3],0])
    
    start_time = time1.time()
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[filters, data_t_for_fit , data_mag_for_fit, data_magerr_for_fit])
    
    nsteps = 200
    
    old_stdout = sys.stdout  
    sys.stdout = open(os.devnull, 'w') # disable printouts
    pos, prob, state = sampler.run_mcmc(starting_positions, nsteps) # לשחק עם זה עד שרואים שהוא מתחיל להתכנס
    sys.stdout = old_stdout # enable printouts
    
    fig, ax = plt.subplots(ndim, sharex=True)
    for i in range(ndim):
        ax[i].plot(sampler.chain[:, :, i].T, '-k', alpha=0.2)
    fig.savefig(chain_filename.replace('chain_','burnin_').replace('.txt','.pdf'))
    
    nsteps = 1000
    
    sampler.reset()
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w') # disable printouts
    pos, prob, state = sampler.run_mcmc(pos, nsteps) 
    sys.stdout = old_stdout # enable printouts
    #print datetime.datetime.now()-n
    fig, ax = plt.subplots(ndim, sharex=True)
    for i in range(ndim):
        ax[i].plot(sampler.chain[:, :, i].T, '-k', alpha=0.2)
    fig.savefig(chain_filename.replace('chain_','walkers_').replace('.txt','.pdf'))
    
    chain_table = Table(sampler.flatchain, names=['M15', 'R500', 'E51', 'toffset_d'])
    chain_table.write(chain_filename,format='ascii.fast_no_header')
    
    samples = sampler.chain[:, 10:, :].reshape((-1, ndim)) #number of steps it throw
    
    results = list(zip(*np.percentile(samples, [16, 50, 84],axis=0)))
    results_line = np.reshape(results,[1,ndim*3])
    print(results_line)
    print('')
    
    s = '{} {} {} {} {} {} {} {} {} {} {} {}'.format(*results_line[0])
    with open(summary_filename,'a') as f: 
        f.write(s)

if loadfromfile:
    
    chain_table = Table.read(filetoload,names=['M15', 'R500', 'E51', 'toffset_d'],format='ascii')

##############################################################################

# fig,axs = beautifulcorner.corner(chain_table)
# fig.savefig(chain_filename.replace('chain_','corner_').replace('.txt','.pdf'))

##############################################################################

offsets = {'UVW2':-3.5, 'UVM2':-2.5, 'UVW1':-1.5, 'U':-0.5, 'B':0.0, 'V':1.8, 'g': 0.8, 'r': 2.5, 'Clear':2.5, 'o': 2.5, 'i': 3.0}

fig,ax = plt.subplots()
t = np.linspace(0,max_days_to_fit*24*3600,50)
for row in random.sample(list(chain_table),k=100):
    old_stdout = sys.stdout  
    sys.stdout = open(os.devnull, 'w') # disable printouts
    mags = calc_model_mags(row['M15'], row['R500'], row['E51'], row['toffset_d'], filters, t)
    sys.stdout = old_stdout # enable printouts
    if mags is not None:
        for f in filters:
            ax.plot(t/(24*3600),np.array(mags[f])+offsets[f],color=filtclr[f],lw=0.5,alpha=0.1)

for f in filters:
    ax.errorbar(data_t[f]/(24*3600),data_mag[f]+offsets[f],yerr=data_magerr[f],fmt='o',color=filtclr[f],mec=filtclr[f],mew=0.5)
    ax.errorbar(data_t_for_fit[f]/(24*3600),data_mag_for_fit[f]+offsets[f],yerr=data_magerr_for_fit[f],fmt='o',color=filtclr[f],mec=filtclr[f],mew=0.5,mfc='w')
    
ax.invert_yaxis()
fig.savefig(chain_filename.replace('chain_','lightcurve_').replace('.txt','.pdf'))

#import corner
#fig = corner.corner(samples[:,:], labels=["$M15$", "$R500$", "$E51$", "$toffset_d$"])
#fig.savefig("sn 2016gkg.png")
#
#
#samples[:, -1] = np.exp(samples[:, -1])
#M15_best, R500_best, E51_best, toffset_d_best = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
#                         zip(*np.percentile(samples, [16, 50, 84],
#                                            axis=0)))
#
#print("M15_best, R500_best, E51_best = ", M15_best, R500_best, E51_best, toffset_d_best)



