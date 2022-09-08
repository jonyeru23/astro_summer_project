# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 10:51:31 2016

@author: iair
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 17:41:48 2016

@author: iair
"""

import os
os.environ['PYSYN_CDBS']='/Users/iair/miniconda2/PYSYN_CDBS/'
import pysynphot
import numpy as np
import astropy
from astropy.table import Table, Column
import emcee, corner
import beautifulcorner
import datetime

from matplotlib import pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
sns.set(style='ticks')
sns.despine    

#figpath = '/Users/iair/Documents/LCOGT/Projects/SN2016gkg/paper/'
figpath = '/Users/iair/Documents/LCOGT/Projects/SN2016gkg/'

##############################################################################

def date2jd(inputdate):

    jd0 = 2451544.5 # On Jan 1, 2000 00:00:00
    td = inputdate-datetime.datetime(2000,1,1,00,00,00)

    return ((td.days*86400)+td.seconds)/float(86400)+jd0

##############################################################################

def get_lcogt_data():
    filepath = '/Users/iair/Documents/LCOGT/Projects/SN2016gkg/'
    filename = 'SN2016gkg_lcophot_1m.txt' # Seems better than the newcal 
    #filename = 'SN2016gkg_lcophot_1m_newcal.txt' 
    data = Table.read(filepath+filename, format='ascii')
    dataout = {}
    dataout['jd'] = data['col1']+2400000.5
    dataout['mag'] = data['col2']
    dataout['mage'] = data['col3']
    dataout['filter'] = data['col4']
    return dataout
    
# Query that gets the data:
#select p.mjd, p.mag, p.dmag, substring(p.filter,1,1) as filter
#from phot as p
#left join photlco as pl on p.photlcoid = pl.id
#where p.targetid = 3218
#and p.issub = 0
#and p.mag < 9999
#and p.dmag < 0.4
#and pl.quality != 1
#order by p.mjd;  

    
##############################################################################

def get_atel_data():
    filepath = '/Users/iair/Documents/LCOGT/Projects/SN2016gkg/'
    filename = 'SN2016gkg_atelphot_9529.txt' 
    dataout = Table.read(filepath+filename, format='ascii', names=['jd','mag','mage','source','filter'])
    return dataout   
    
##############################################################################

def get_swift_data():
    filepath = '/Users/iair/Documents/LCOGT/Projects/SN2016gkg/'
    filename = 'SN2016gkg_swiftphot.txt' 
    data = Table.read(filepath+filename, format='ascii', names=['date','jd','B','Be','UVM2','UVM2e','UVW1','UVW1e','U','Ue','UVW2','UVW2e','V','Ve','one','tel','filt'])
    dataout = {'jd':data['jd'],'mag':np.array([]),'mage':np.array([]),'filter':np.array([])}
    for i,row in enumerate(data):
        dataout['mag'] = np.append(dataout['mag'],data[row['filt']][i])
        dataout['mage'] = np.append(dataout['mage'],data[row['filt']+'e'][i])
        dataout['filter'] = np.append(dataout['filter'],row['filt'])
    return dataout  

##############################################################################

def get_swift_data_pb():
    filepath = '/Users/iair/Documents/LCOGT/Projects/SN2016gkg/'
    filename = 'SN2016gkg_swiftphot_peterbrown.txt' 
    dataout = Table.read(filepath+filename, format='ascii', names=['filter','jd','mag','mage','3SigMagLim', '0.98SatLim', 'Rate', 'RateErr','Aperture','Frametime','Exp','Telapse'])
    dataout['jd'] = dataout['jd']+2400000.5
    return dataout

##############################################################################

def get_sv_data():
    filepath = '/Users/iair/Documents/LCOGT/Projects/SN2016gkg/'
    filename = 'SN2016gkg_svphot.txt' 
    skip = True
    filts = []
    jds = []
    mags = []
    dmags = []
    nondets = []
    instrs = []
    with open(filepath+filename) as f:
        for line in f:
            if line[0]=='#': continue
            if line[0]=='_': skip = True
            vals = line.split()
            if not vals: continue
            if vals[0]=='Date':
                skip = False
                filters = vals[2:-2]
                continue
            if not skip:
                jd = float(vals[1])
                instr = vals[-1]
                src = int(vals[2*len(filters)+2])
                nondet = src < 0
                if nondet: err = np.nan
                for i,filt in enumerate(filters):
                    mag = float(vals[2*i+2])
                    if not nondet: err = float(vals[2*i+3])
                    if mag < 9999:
                        filts.append(filt)
                        jds.append(jd)
                        mags.append(mag)
                        dmags.append(err)
                        nondets.append(nondet)
                        instrs.append(instr)
    datatable = Table([jds, filts, mags, dmags, nondets, instrs], names=['jd', 'filter', 'mag', 'mage', 'nondet', 'instr'])
    return datatable
    
##############################################################################

def calc_extinction(Av,Rv):
        x = np.array([0.1928, 0.2246, 0.26, 0.34, 0.43, 0.54, 0.64, 0.80, 0.36, 0.47, 0.62, 0.75, 0.89, 1.25, 1.66, 2.19])**-1 # wavenumber in inverse microns
        y = x - 1.82
        aO = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
        bO = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
        aIR = 0.574*x**1.61
        bIR = -0.527*x**1.61
        aUV = 1.752 - 0.316*x - 0.104/((x - 4.67)**2 + 0.341)  # only valid above 169 nm
        bUV = -3.090 + 1.825*x + 1.206/((x - 4.62)**2 + 0.263) # only valid above 169 nm
        a = aO.copy()
        b = bO.copy()
        a[x < 1.1] = aIR[x < 1.1]
        b[x < 1.1] = bIR[x < 1.1]
        a[x > 3.3] = aUV[x > 3.3]
        b[x > 3.3] = bUV[x > 3.3]
        A = Av * (a + b / Rv)
        extinction = {f: Af for f, Af in zip('2M1UBVRIugrizJHK', A)}
        extinction['0'] = extinction['R']
        # 2 = UVW2, M = UVM2, 1 = UVW1, and 0 = unfiltered
        return extinction
        
##############################################################################

def correct_for_extinction(data,extinction):
    
    for i,row in enumerate(data['mag']):
        data['mag'][i] = data['mag'][i]-extinction[data['filter'][i]]
    return data

##############################################################################

def bin_daily(data):
    
    bdata = {'jd':np.array([]),'mag':np.array([]),'mage':np.array([]),'filter':np.array([]),'source':np.array([])}   
    
    filters = set(data['filter'])
    data = Table(data)
    if 'source' not in data.colnames:
        data['source'] = np.tile('s',len(data['jd']))
    sources = set(data['source'])    
    
    for filt in filters:
        for source in sources:
        
            w = (data['filter']==filt) & (data['source']==source)
            if sum(w) > 0:
        
                jd = data['jd'][w]
                mag = data['mag'][w]
                mage = data['mage'][w]
                                
                # Sort:
                s = np.argsort(jd)
                jd = jd[s]
                mag = mag[s]
                mage = mage[s]
        
                # Convert to flux:
                flux = 10**(-0.4*mag)
                fluxe = abs(flux*(-0.4)*np.log(10)*mage)
                if 0 in fluxe:
                    weight = np.tile(1,len(fluxe))
                else:
                    weight = 1/fluxe**2
                
                # Do the binning using weighted averages:
                newdayjump = 0.5
                toaverage = [0]
                binned_jd = np.array([])
                binned_flux = np.array([])
                binned_fluxe = np.array([])
                for i in range(0,len(jd)):
                    if abs(jd[i]-jd[toaverage[-1]]) < newdayjump:
                        toaverage.append(i)
                    else:
                        binned_jd = np.append(binned_jd,sum(weight[toaverage]*jd[toaverage])/sum(weight[toaverage]))
                        binned_flux = np.append(binned_flux,sum(weight[toaverage]*flux[toaverage])/sum(weight[toaverage]))
                        binned_fluxe = np.append(binned_fluxe,np.sqrt(1/sum(weight[toaverage])))
                        toaverage = [i]
                    if i == len(jd)-1: # we've reached the end
                        binned_jd = np.append(binned_jd,sum(weight[toaverage]*jd[toaverage])/sum(weight[toaverage]))
                        binned_flux = np.append(binned_flux,sum(weight[toaverage]*flux[toaverage])/sum(weight[toaverage]))
                        binned_fluxe = np.append(binned_fluxe,np.sqrt(1/sum(weight[toaverage])))
                    
                # Convert back to magnitudes:
                binned_mag = -2.5*np.log10(binned_flux)
                if 0 in fluxe:
                    binned_mage = np.tile(0,len(binned_flux))
                else:
                    binned_mage = abs(-2.5*binned_fluxe/(binned_flux*np.log(10)))
                    
                bdata['jd'] = np.concatenate((bdata['jd'],binned_jd))
                bdata['mag'] = np.concatenate((bdata['mag'],binned_mag))
                bdata['mage'] = np.concatenate((bdata['mage'],binned_mage))
                bdata['filter'] = np.concatenate((bdata['filter'],np.tile(filt,len(binned_jd))))
                bdata['source'] = np.concatenate((bdata['source'],np.tile(source,len(binned_jd))))
        
    s = np.argsort(bdata['jd'])
    bdata['jd'] = bdata['jd'][s]
    bdata['mag'] = bdata['mag'][s]
    bdata['mage'] = bdata['mage'][s]
    bdata['filter'] = bdata['filter'][s]
    bdata['source'] = bdata['source'][s]
        
    return Table(bdata)

##############################################################################

def get_filtclr():

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
    
    return filtclr

##############################################################################

def export_data(data,filename):
    
    colnames = ['jd','filter','mag','mage']
    t = Table([data[k] for k in colnames],names=colnames)
    t.write(filename,format='ascii')

##############################################################################

def plot_lc(re=None,me=None,ve=None,toffset=None,bersten=False):
      
    logplot = False
    colors = get_filtclr() 
    fs = 14
    ms = 5
    lw = 1
    #yul=16.7
    #yll=22

    fig,ax = plt.subplots(figsize=(10,8))
    #fig,ax = plt.subplots()    
    
    plotdata = {'jd':[],'telescope':[],'source':[],'filter':[],'mag':[],'mage':[],'islim':[]}    
    
    # MODELS

    if re is not None and ve is not None and me is not None and toffset is not None:       
       
        lss = ['--','-.','-']
        lws = [1,1,0.5]
        for j in range(len(re)):
        
            if loadfromfile != '':
                thismodel = loadfromfile[j].split('_')[0]
                thisn = float(loadfromfile[j].split('_')[3].split('e')[0])
            else:
                thismodel = model
                thisn = n
            if j==0:
                tstart = 0.001
            else:
                tstart = 0.001
            t = np.linspace(tstart,tend_fit-toffset[j],200)*86400
            if thismodel == 'piro':
                lbol, radius, temperature = eval_piro(re[j],me[j],ve[j],t)
            if thismodel =='sapir':
                lbol, radius, temperature, tdmin, tdmax = eval_sapir(re[j],me[j],ve[j],thisn,t)
                w = (t-toffset[j]>=tdmin*86000) & (t-toffset[j]<=tdmax*86000)
                t = t[w]
                lbol = lbol[w]
                radius = radius[w]
                temperature = temperature[w]
                #print tdmin,tdmax
            radius_in_rsun = radius / astropy.constants.R_sun.to('cm').value

            if fit_ni:
                ni_fits = fit_nickel(lcogt1m_data,['B','V','g','r','i'],toffset[j])
                ni_fits_swift = fit_nickel(swift_data,['U','UVW1','UVW2','UVM2'],toffset[j])
                for filt in ['U','UVW1','UVW2','UVM2']:
                    w = ni_fits_swift['filter'].index(filt)
                    ni_fits['filter'].append(filt)
                    ni_fits['polynomial'].append(ni_fits_swift['polynomial'][w])
    
            for filtname in ['swift,UVW2','swift,UVM2','swift,UVW1','U','B','V','g','r','i','o']:

                if 'swift' in filtname:
                    clrfiltname = filtname.split('swift,')[1]
                else:
                    clrfiltname = filtname                
                
                modelmag = []
                if filtname in ['g','r','i']: 
                    newfiltname = 'sdss,'+filtname
                else:
                    newfiltname = filtname 
                bp = bandpasses[newfiltname]                 
                                
                if fit_ni:
                    ni_flux = None
                    if clrfiltname in ni_fits['filter']:
                        w = ni_fits['filter'].index(clrfiltname)
                        if ni_fits['polynomial'][w] is not None:
                            ni_flux = np.polyval(ni_fits['polynomial'][w],t/86400.0-toffset[j])                                        
                
                for i,x in enumerate(t):
                    bb = pysynphot.BlackBody(temperature[i]) # For 1 Solar Radius at 1 kpc 
                    obs = pysynphot.Observation(bb, bp)
                    if 'sdss' in newfiltname:
                        magsystem = 'ABMag'
                    else:
                        magsystem = 'VegaMag'
                    expected_mag = obs.effstim(magsystem) - 2.5*np.log10((radius_in_rsun[i]**2)*((1000.0/distance)**2)) # Rescaling from the default (1 solar radius at 1000 pc)
                    if fit_ni:                    
                        if ni_flux is not None:
                            expected_flux = 10**(-0.4*expected_mag)+ni_flux[i]
                            expected_mag = -2.5*np.log10(expected_flux)
                    modelmag.append(expected_mag)
                    #if expected_mag > 18:
                    #    break                
                if clrfiltname in offsets:
                    offset = offsets[clrfiltname]
                else:
                    offset = 0
                ax.plot(np.array(t[0:len(modelmag)]/86400.0)+toffset[j],np.array(modelmag)+offset,color=colors[clrfiltname],lw=lws[j], ls=lss[j])
#                if ni_flux is not None:
#                    ax.plot(np.array(t[0:len(t)]/86400.0)+toffset[j],-2.5*np.log10(np.array(ni_flux))+offset,color=colors[clrfiltname],lw=lws[j], ls=lss[j])                    
        
    if bersten:
        
        bmodel = Table.read('bersten_2011dh_models_r270.csv',format='ascii',names=['day','gmag'])
        ax.plot(bmodel['day']-toffset,bmodel['gmag']+dm+offsets['g'],color=colors['g'],lw=lw+1,ls='--')

    # DATA

    #PROMPT Data:    
    data = dlt40_data
    for filt in ['Clear']:
        w = (data['filter']==filt) & ((data['jd']-t0)/(1+redshift)<tend_plot)
        mec = 'k'
        mew = 0.5
        #ax.errorbar((data['jd'][w]-t0)/(1+redshift)-toffset,data['mag'][w],yerr=data['mage'][w],fmt='p',color=mec,mfc=colors[filt],mec=mec,mew=mew,label='Prompt '+filt,ms=ms+2)    

    #LCOGT Data:
    data = lcogt1m_data
    
    # Simulate the error weighting from the fit:
#    t = (data['jd']-t0)/(1+redshift)-toffset
#    flux = 10**(-0.4*data['mag'])
#    fluxe = abs(flux*(-0.4)*np.log(10)*data['mage'])
#    minfluxe = abs(flux*(-0.4)*np.log(10)*0.01)
#    newfluxe = [max([fluxe[i]*t[i]/tend_fit,minfluxe[i]]) for i in range(len(fluxe))]
#    weightedmage = abs(newfluxe/(-0.4*np.log(10)*flux))
    
    for filt in ['B','V','g','r','i']:
        w = (data['filter']==filt) & ((data['jd']-t0)/(1+redshift)<tend_plot)
        if filt in offsets:
            offset = offsets[filt]
        else:
            offset = 0
        ax.errorbar((data['jd'][w]-t0)/(1+redshift),data['mag'][w]+offset,yerr=data['mage'][w],fmt='o',color=colors[filt],mec=colors[filt],mew=0.5,label='LCOGT-1m '+filt,ms=ms)
        for key in ['jd','mag','mage']:       
            plotdata[key].extend(list(data[key][w]))
        plotdata['filter'].extend(list(np.tile(filt,np.sum(w))))
        plotdata['islim'].extend(list(np.tile(False,np.sum(w))))
        plotdata['telescope'].extend(list(np.tile('LCOGT',np.sum(w))))
        plotdata['source'].extend(list(np.tile('This Work',np.sum(w))))
        
#        w = (data['filter'] == filt) & (data['jd']-t0 >= tstart_ni) & (data['jd']-t0 <= tend_plot)
#        x = data['jd'][w]-t0-toffset
#        x = np.append(0,x)
#        y = 10**(-0.4*data['mag'][w])
#        y = np.append(0,y)
#        p = np.polyfit(x,y,2)
#        xx = np.linspace(0,tend_plot,50)
#        yy = np.polyval(p,xx)
#        ax.plot(xx,-2.5*np.log10(yy)+offset,ls='-',color=colors[filt])
    
    #ATel Data:
    fmts = {'Buso&Otero':'*','ASAS-SN':(6,0,90),'ATLAS':(5,0,0),'LCOGT':'s'}
    outnames = {'Buso&Otero':'Buso & Otero','ASAS-SN':'ASAS-SN','ATLAS':'ATLAS','LCOGT':'LCOGT'}
    outsources = {'Buso&Otero':'TNS','ASAS-SN':'Chen et al. 2016','ATLAS':'This Work','LCOGT':'Chen et al. 2016'}   
    data = atel_data
    
    # Simulate the error weighting from the fit:
#    t = (data['jd']-t0)/(1+redshift)-toffset
#    flux = 10**(-0.4*data['mag'])
#    fluxe = abs(flux*(-0.4)*np.log(10)*data['mage'])
#    minfluxe = abs(flux*(-0.4)*np.log(10)*0.01)
#    newfluxe = [max([fluxe[i]*t[i]/tend_fit,minfluxe[i]]) for i in range(len(fluxe))]
#    weightedmage = abs(newfluxe/(-0.4*np.log(10)*flux))

    for filt in ['Clear','B','V','g','r','i','o']:
        for src in fmts:
            if filt in offsets:
                offset = offsets[filt]
            else:
                offset = 0
            w = (data['filter']==filt) & ((data['jd']-t0)/(1+redshift)<tend_plot) & (data['source']==src)
            if filt == 'Clear':
                mec = 'k'
                mew = 0.5
            else:
                mec = colors[filt]
                mew = 0
            if fmts[src] in [(6,0,90),(5,0,0)]:
                thisms = ms+2
            elif fmts[src] == '*':
                thisms = ms+3
            else:
                thisms = ms
            ax.errorbar((data['jd'][w]-t0)/(1+redshift),data['mag'][w]+offset,yerr=data['mage'][w],ls='none',marker=fmts[src],color=mec,mfc=colors[filt],mec=mec,mew=0.5,label='ATel '+filt,ms=thisms)
            for key in ['jd','mag','mage']:       
                plotdata[key].extend(list(data[key][w]))
            plotdata['filter'].extend(list(np.tile(filt,np.sum(w))))
            plotdata['islim'].extend(list(np.tile(False,np.sum(w))))
            plotdata['telescope'].extend(list(np.tile(outnames[src],np.sum(w))))
            plotdata['source'].extend(list(np.tile(outsources[src],np.sum(w))))
    
    
    #PESSTO Data:    
    data = pessto_data
    for filt in ['U','B','V','R','I']:
        w = (data['filter']==filt) & ((data['jd']-t0)/(1+redshift)<tend_plot)
        if filt in offsets:
            offset = offsets[filt]
        else:
            offset = 0
        #ax.errorbar((data['jd'][w]-t0)/(1+redshift),data['mag'][w],yerr=data['mage'][w],fmt='h',color=colors[filt],label='EFOSC '+filt,ms=ms+2)    
        
    #Swift Data:
    data = swift_data
    for filt in ['UVW1','UVM2','UVW2','U','B','V']:
        if filt in offsets:
            offset = offsets[filt]
        else:
            offset = 0
        w = (data['filter']==filt) & ((data['jd']-t0)/(1+redshift)<tend_fit)
        if filt in offsets:
            offset = offsets[filt]
        else:
            offset = 0
        ax.errorbar((data['jd'][w]-t0)/(1+redshift),data['mag'][w]+offset,yerr=data['mage'][w],marker='D',color=colors[filt],ls='none',mfc=colors[filt],mec=colors[filt],mew=0.5,label='Swift '+filt,ms=ms)
        for key in ['jd','mag','mage']:       
            plotdata[key].extend(list(data[key][w]))
        plotdata['filter'].extend(list(np.tile(filt,np.sum(w))))
        plotdata['islim'].extend(list(np.tile(False,np.sum(w))))
        plotdata['telescope'].extend(list(np.tile('Swift',np.sum(w))))
        plotdata['source'].extend(list(np.tile('This Work',np.sum(w))))
    
    #Amateurs:
    #ax.plot((2457652.79-t0)/(1+redshift)-toffset,15+offsets['V'],marker='*',ls='none',color=colors['V'])
    
    #Nicholls:
    data = nicholls_data
    filt = nicholls_data['filter'][0]
    if filt in offsets:
        offset = offsets[filt]
    else:
        offset = 0
    ax.errorbar((data['jd'][0]-t0)/(1+redshift),data['mag'][0]+offset,yerr=data['mage'][0],marker=(6,1,0),mfc=colors[filt],label='Nicholls '+filt,ms=ms+2,mec='k',color='k',mew=0.5)
    for key in ['jd','mag','mage']:       
        plotdata[key].append(data[key][0])
    plotdata['filter'].append(filt)
    plotdata['islim'].append(False)
    plotdata['telescope'].append('Nicholls')
    plotdata['source'].append('This Work')
    
    #Monard:
    data = monard_data
    filt = monard_data['filter'][0]
    if filt in offsets:
        offset = offsets[filt]
    else:
        offset = 0
    ax.errorbar((data['jd'][0]-t0)/(1+redshift),data['mag'][0]+offset,yerr=data['mage'][0],marker=(7,1,0),mfc=colors[filt],label='Monard '+filt,ms=ms+2,mec='k',color='k',mew=0.5)    
    for key in ['jd','mag','mage']:       
        plotdata[key].append(data[key][0])
    plotdata['filter'].append(filt)
    plotdata['islim'].append(False)
    plotdata['telescope'].append('Monard')
    plotdata['source'].append('This Work')       
    
    #Sanchez:
    data = sanchez_data    
    filt = nicholls_data['filter'][0]
    if filt in offsets:
        offset = offsets[filt]
    else:
        offset = 0
    ax.errorbar((data['jd'][0]-t0)/(1+redshift),data['mag'][0]+offset,yerr=data['mage'][0],marker=(5,1,180),mfc=colors[filt],label='Sanchez '+filt,ms=ms+2,mec='k',color='k',mew=0.5)
    for key in ['jd','mag','mage']:       
        plotdata[key].append(data[key][0])
    plotdata['filter'].append(filt)
    plotdata['islim'].append(False)
    plotdata['telescope'].append('Sanchez')
    plotdata['source'].append('This Work')

    #ax.legend(ncol=2, numpoints=1, frameon=False, loc='best',fontsize=fs-2)      
    
    ax.invert_yaxis()
    ax.tick_params(axis='x', labelsize=fs)
    ax.tick_params(axis='y', labelsize=fs)
    ax.xaxis.set_tick_params(width=2, length=7)
    ax.yaxis.set_tick_params(width=2, length=7)
    if logplot:
        ax.set_xlim([0.2,tend_plot+2])
        ax.set_xscale('log')
        ax.set_ylim([21,9.7])
    else:
        ax.set_xlim([-0.5,tend_plot])
        ax.set_ylim([21,8])    
    ax.set_xlabel('Time since discovery (rest-frame days)', fontsize=fs)
    #ax.set_xlabel('Time since explosion (rest-frame days)', fontsize=fs)
    [i.set_linewidth(1.5) for i in ax.spines.itervalues()]
    plt.ylabel('Apparent magnitude', fontsize=fs)
    
    #ATel Limit (need ot plot after reversing axis):
    filt = 'V'
    x = 2457651.6653
    y = 17.36-extinction[filt]
    ax.errorbar((x-t0)/(1+redshift),y+offsets[filt],yerr=0.5,lolims=True,ls='',capsize=5,color=colors[filt])
    ax.plot([(x-t0)/(1+redshift)-0.05,(x-t0)/(1+redshift)+0.05],[y+offsets[filt],y+offsets[filt]],color=colors[filt])
    
    plotdata['jd'].append(x)
    plotdata['mag'].append(y)
    plotdata['mage'].append(0)
    plotdata['filter'].append(filt)
    plotdata['islim'].append(True)
    plotdata['telescope'].append('ASAS-SN')
    plotdata['source'].append('Nicholls et al. 2016')    
    
    plottable = Table([plotdata[field] for field in ['jd','telescope','source','filter','mag','mage','islim']], names=['JD','Telescope','Source','Filter','Magnitude','Magnitude_Error','Is_Limit'])
    plottable.sort('JD')
    plottable.columns['JD'].format='0.3f'
    plottable.columns['Magnitude'].format='0.3f'
    plottable.columns['Magnitude_Error'].format='0.3f'
    print(plottable)
    plottable.write('fig1data_extinction_corrected_final.txt',format='ascii')

    # Abs Mag Axis
    ax2 = plt.twinx()
    ax2.set_ylim([m-dm for m in ax.get_ylim()])
    ax2.tick_params(axis='y', labelsize=fs)
    ax2.yaxis.set_tick_params(width=2, length=7)
    [i.set_linewidth(1.5) for i in ax2.spines.itervalues()]
    ax2.set_ylabel('Absolute magnitude', fontsize=fs)
     
    #Legend:
    #telescopes = ['Buso & Otero','ATLAS','ASAS-SN','LCOGT (ATel)','LCOGT','NTT','PROMPT','Swift']
    #markers = ['*',(5,0,180),(6,0,90),'s','o','h','p','D']
    telescopes = ['Buso','Sanchez','Nicholls','Monard','ASAS-SN','LCOGT (C16)','LCOGT','ATLAS','Swift']
    markers = ['*',(5,1,180),(6,1,0),(7,1,0),(6,0,90),'s','o',(5,0,0),'D']

    telescope_lines = []
    for i,telescope in enumerate(telescopes):
        if markers[i] == '*':
            corr = 3
        elif markers[i] == 'p' or markers[i] == 'h' or markers[i] == (5,0,0) or markers[i] == (6,0,90) or markers[i] == (6,1,0) or markers[i] == (7,1,0) or markers[i] == (5,1,180):
            corr = 2
        else:
            corr = 0
        if '-' in markers[i]:
            m = plt.Line2D([], [], linestyle=markers[i], color='k', ms=ms+1+corr)
        else:
            m = plt.Line2D([], [], linestyle='none', marker=markers[i], mec='k', mfc='w', mew=1, ms=ms+1+corr)
        telescope_lines.append(m)
    if logplot:
        ax.legend(telescope_lines, telescopes, ncol=2,numpoints=1, frameon=False, fontsize=fs-2, loc='lower right')
    else:
        ax.legend(telescope_lines, telescopes, ncol=2,numpoints=1, columnspacing=1.5, frameon=False, fontsize=fs-2, loc=[0.59,0.62])
   
    #filters = ['g','r','i','o','U','UVW1','UVM2','UVW2','B','V','R','I','Clear']
    #filterlabels = ['g','r','i','o','U','UVW1','UVM2','UVW2','B','V','R','I','Clear']
    filters = ['UVW1','UVM2','UVW2','U','g','r','i','o','B','V','Clear']
    filterlabels = ['UVW1','UVM2','UVW2','U','g','r','i','o','B','V','Clear']
    for i,filt in enumerate(filterlabels):
        if filt in offsets:
            if offsets[filt] != 0:
                offsetstr = '%+0.1f' %offsets[filt]
                filterlabels[i] += offsetstr.replace('-',' - ').replace('+',' + ')
    
    filt_lines = []
    for filt in filters:
        if filt == 'Clear':
            m = Patch(ec='k', fc='w', lw=1)
        else:
            m = Patch(ec=colors[filt], fc=colors[filt], lw=1)
        filt_lines.append(m) 
    if logplot:                   
        ax2.legend(filt_lines, filterlabels, ncol=2, numpoints=1, frameon=False, loc='upper right', fontsize=fs-2)
        fnpostfix = '_log'
    else:
        ax2.legend(filt_lines, filterlabels, ncol=3, numpoints=1, columnspacing=1.5, frameon=False, loc=[0.44,0.82], fontsize=fs-2)
        fnpostfix = '_linear'
    if fit_ni:
        fnpostfix = '_ni'
    else:
        fnpostfix = ''

    if re:
        ax3 = plt.twiny()
        ax3.tick_params(axis='x', labelsize=0)
        ax3.xaxis.set_tick_params(width=0, length=0)
        models = ['P15','SW16 (n = 3/2)','SW16 (n = 3)']
        markers = lss

        model_lines = []
        for i,mdl in enumerate(models):
            m = plt.Line2D([], [], linestyle=markers[i], color='k')
            model_lines.append(m)
        ax3.legend(model_lines, models, ncol=1,frameon=False, fontsize=fs-2, loc=[0.765,0.49])    

    if re:
        fig.savefig(figpath+'photometry_'+model+'_'+weighting+'_'+velocity+'_'+str(n)+fnpostfix+'.pdf', bbox_extra_artists=(ax2,), bbox_inches='tight')
    else:
        fig.savefig(figpath+'photometry_'+weighting+'_'+velocity+'_'+str(n)+fnpostfix+'.pdf', bbox_extra_artists=(ax2,), bbox_inches='tight')

##############################################################################

def eval_piro(re,me,ve,t):

    lbol = 8.27e11*re*(ve**2)*np.exp((-4.135e-20)*t*(t+2*re/ve)*ve*((me/0.01)**(-1)))

    radius = re+ve*t
    temperature = (lbol / (4*np.pi*radius**2*sigma))**(0.25)
    
    #print lbol,radius,temperature    
    
    return (lbol,radius,temperature)
    
##############################################################################

def eval_sapir(r,me,v,n,t,mc=1):    
    
    m = mc+me # Solar masses
    td = t/86400.0 # t in sec, td in days  
    r13 = r/1.0e13
    v85 = v/(10**8.5)
    
    if n == 3/2.0:        
        f = (me/mc)**0.5 
        tpref = 1.61
        tcolpref = 1.1
        lpref = 2
        eps1 = 0.027
        eps2 = 0.086
        capitala = 0.94
        a = 1.67
        alpha = 0.8
    if n == 3:
        f = 0.08*me/mc
        tpref = 1.69
        tcolpref = 1.0
        lpref = 2.1
        eps1 = 0.016
        eps2 = 0.175
        capitala = 0.79
        a = 4.57
        alpha = 0.73
        
    tdmin = 0.2*r13/v85*np.max([0.5,r13**0.4/((f*m)**0.2*v85**0.7)])
    tdmax = 7.4*r13**0.55
    ttrd = 19.5*(me/v85)**0.5
    
    lrw = lpref * 1e42 * (v85*td**2/(f*m))**(-eps2) * v85**2*r13
    l = lrw * capitala * np.exp(-(a*td/ttrd)**alpha)
    
    Tph = tpref * ((v85**2)*(td**2)/(f*m))**eps1 * r13**0.25 * td**(-0.5) / kb
    Tcol = tcolpref * Tph
        
    lbol = l
    temperature = Tcol
    radius = (lbol / (4*np.pi*sigma*temperature**4))**0.5
    
    return (lbol,radius,temperature,tdmin,tdmax)

##############################################################################

def fit_nickel(data,filternames,toffset):
    
    ni_fits = {'filter':[],'polynomial':[]}
    for filtername in filternames:
        ni_fits['filter'].append(filtername)
        w = (data['filter'] == filtername) & (data['jd']-t0 >= tstart_ni) & (data['jd']-t0 <= tend_plot)
        x = data['jd'][w]-t0-toffset
        x = np.append(0,x)
        y = 10**(-0.4*data['mag'][w])
        y = np.append(0,y)
        if len(x)>3:
            p = np.polyfit(x,y,2)
            ni_fits['polynomial'].append(p)
        else:
            ni_fits['polynomial'].append(None)  
    return ni_fits

##############################################################################

def log_flat_prior(theta):
    re13, ve9, me, toffset = theta
    # re13 in 10^13 cm
    # ve9 in 10^9 cm/s
    # me in solar masses
    # toffset in days relative to discovery (t0)
    re = re13 * 1e13
    ve = ve9 * 1e9
    if re > 1e10 and re < 1e14 and ve > 100e5 and ve < 100000e5 and me > 0.005 and me < 1 and toffset > -0.5 and toffset < 0:
        return 0 # log(1)
    else:
        return -np.inf  # log(0)
    
##############################################################################    
    
def log_logflat_prior(theta):
    if velocity == 'fixed':
        re13, me, toffset = theta
    else:
        re13, me, ve9, toffset = theta
    # re13 in 10^13 cm
    # ve9 in 10^9 cm/s
    # me in solar masses
    # toffset in days relative to discovery (t0)
    re = re13 * 1e13
    if velocity == 'free':
        ve = ve9 * 1e9
        if re > 1e10 and re < 1e14 and ve > 100e5 and ve < 100000e5 and me > 0.001 and me < 1 and toffset > -0.5 and toffset < 0:
            #return (-np.log10(re)) + (-np.log10(ve)) + (-np.log10(me)) + (-np.log10(-toffset))
            return (-np.log10(re)) + (-np.log10(ve)) + (-np.log10(me))
        else:
            return -np.inf  # log(0)
    else:
        if re > 1e10 and re < 1e14 and me > 0.001 and me < 1 and toffset > -0.5 and toffset < 0:
            #return (-np.log10(re)) + (-np.log10(ve)) + (-np.log10(me)) + (-np.log10(-toffset))
            return (-np.log10(re)) + (-np.log10(me))
        else:
            return -np.inf  # log(0)

##############################################################################

def log_likelihood(theta, filternames, bandpasses, all_t, all_m, all_merr, distance):
    # m in apparent mag

    if velocity == 'fixed':
        re13, me, toffset = theta
    else:
        re13, me, ve9, toffset = theta
    # re13 in 10^13 cm
    # ve9 in 10^9 cm/s
    # me in solar masses
    # toffset in days relative to discovery (t0) 
    re = re13 * 1e13

    if velocity == 'free':
        ve = ve9 * 1e9        
        if re < 1e10 or re > 1e14 or ve < 100e5 or ve > 100000e5 or me < 0.001 or me > 1 or toffset < -0.5 or toffset > 0:
            return -np.inf    
    else:
        ve = fixedve
        if re < 1e10 or re > 1e14 or me < 0.001 or me > 1 or toffset < -0.5 or toffset > 0:
            return -np.inf    
            
    #plot_lc(re,ve,me,toffset)
            
    # Fit the Ni curve:
    if fit_ni:
        ni_fits = fit_nickel(lcogt1m_data,filternames,toffset)
        for filtername in filternames:
            w = ni_fits['filter'].index(filtername)
            if ni_fits['polynomial'][w] is None: # No data form LCO, try Swift
                ni_fits_swift = fit_nickel(swift_data,[filtername],toffset)
                if ni_fits_swift['polynomial'][0] is not None: # Found in Swift!
                    ni_fits['polynomial'][w] = ni_fits_swift['polynomial'][0]
    
    chi2 = []
    normalization = []
    expected_mags = []
    
    for filtername in filternames:
        
        t = all_t[filtername]
        m = all_m[filtername]
        merr = all_merr[filtername]

        if model == 'piro':
            lbol, radius, temperature = eval_piro(re,me,ve,(t-toffset)*86400)
        if model == 'sapir':
            lbol, radius, temperature, tdmin, tdmax = eval_sapir(re,me,ve,n,(t-toffset)*86400)
        radius_in_rsun = radius / astropy.constants.R_sun.to('cm').value
        
        # Ni contribution:
        if fit_ni:
            w = ni_fits['filter'].index(filtername)
            if ni_fits['polynomial'][w] is not None:
                ni_flux = np.polyval(ni_fits['polynomial'][w],t-toffset)
            else:
                ni_flux = None

        bp = bandpasses[filtername]   

        for i,x in enumerate(t):

            if (t[i]==0 and model == 'sapir') or np.isnan(temperature[i]):
                continue
            bb = pysynphot.BlackBody(temperature[i]) # For 1 Solar Radius at 1 kpc 
            obs = pysynphot.Observation(bb, bp)
            if ('sdss' in filtername) or filtername == 'o':
                magsystem = 'ABMag'
            else:
                magsystem = 'VegaMag'
            try:
                expected_mag = obs.effstim(magsystem) - 2.5*np.log10((radius_in_rsun[i]**2)*((1000.0/distance)**2)) # Rescaling from the default (1 solar radius at 1000 pc)
            except:
                expected_mag = 50
            expected_mags.append(expected_mag)

            # Convert to flux for chi2 measurement:
            measured_flux = 10**(-0.4*m[i])
            measured_fluxe = abs(measured_flux*(-0.4)*np.log(10)*merr[i])
            if weighting == 'early':
                # Weigh by time to give more weight (eff. smaller error) to early data - but don't lower error below 0.005 mags:
                minmage = 0.005
                minfluxe = abs(measured_flux*(-0.4)*np.log(10)*minmage)
                measured_fluxe = max([measured_fluxe*x/tend_fit,minfluxe])
            expected_flux = 10**(-0.4*expected_mag)
            if fit_ni:            
                if ni_flux is not None:
                    expected_flux = expected_flux + ni_flux[i]
            
            thischi2 = ((expected_flux-measured_flux)/measured_fluxe)**2
            chi2.append(thischi2)
            normalization.append(np.log(2 * np.pi * measured_fluxe**2))

    #print chi2, normalization
    #print -0.5 * np.sum(normalization + chi2)
    return -0.5 * np.sum(normalization + chi2)
    
##############################################################################

def log_posterior(theta, filternames, bandpasses, all_t, all_m, all_merr, distance):
    return log_logflat_prior(theta) + log_likelihood(theta, filternames, bandpasses, all_t, all_m, all_merr, distance)
    
##############################################################################

# Piro Ni Early Free v looks good but parameter distributions are multi-peaked...
# Much better without the two new amateur points without Ni, set v to 1.7, same for Sapir with n=3/2.

# Sapir 3/2 weighted, no Ni, fixed to v=1.7 looks ok. Try unweighted.
# Sapir 3/2 weighted, no Ni, velocity free looks ~ok. Try with more steps. Tried with 4000, fit good, dstirbution meh.

models = ['piro','sapir']
model = models[1]

fit_ni = False

ns = [3/2.0,3] # n=3<->Radiative,BSG; n=3/2<->Convective,RSG
n = ns[1]

weightings = ['normal','early']
weighting = weightings[1] 

velocities = ['free','fixed']
velocity = velocities[0]

fixedve9 = 1.7 #~17000 km/s from SALT spectrum reported in ATel
fixedve = fixedve9 * 1e9

nsteps = 4000

# New Tartaglia+ host extinction estimates:
# E(B-V) D2 =  0.0993573500587  [ 0.0703396171216 0.146960289898 ]
# E(B-V) D1 =  0.0855460576341  [ 0.0578362331783 0.132495155746 ]

ebvs = [0.09,0.06,0.15] #avg D1 & D2, lower D1, upper D2 
ebv = ebvs[1]

#t0 = 2457651.7484 # TNS Discovery JD (average of first images).
t0 = 2457651.69028 # VSNET Discovery JD (actual first image).
tend_fit = 3.3
tend_plot = 6
tstart_ni = 5 # For removing the Ni contribution
redshift = 0.00494
#dm = 31.97 # NED
#distanc)e = 25.267e6 # NED, in pc
dm = 32.11 # Tully-Fisher measurements \citep{Tully09}.
distance = 26.4e6 # Tully-Fisher measurements \citep{Tully09}, in pc

kb = astropy.constants.k_B.to('eV/K').value
sigma = astropy.constants.sigma_sb.to('erg s^-1 cm^-2 K^-4').value   

#offsets = {'B':0.4, 'V':-0.5, 'g': 0.2, 'r': -0.6, 'i': -0.2, 'UVW1':-0.2, 'UVM2':-1, 'UVW2':-0.8}
#offsets = {key:offsets[key]*3 for key in offsets}
offsets = {'UVW2':-3.5, 'UVM2':-2.5, 'UVW1':-1.5, 'U':-0.5, 'B':0.0, 'V':1.8, 'g': 0.8, 'r': 2.5, 'Clear':2.5, 'o': 2.5, 'i': 3.0}

# Extinction:

r_v = 3.1

# MW Extinction

mw_extinction = {'U':0.084, 'B':0.070, 'V':0.053, 'R':0.042, 'I':0.029, 'u':0.082, 'g':0.064, 'r':0.044, 'i':0.033, 'z':0.025} # Schlafly+ 2011 via NED
mw_extinction['Clear'] = mw_extinction['R']
mw_extinction['o'] = mw_extinction['R']

# Calculate for Swift UV Filters
calculated_extinctions = calc_extinction(mw_extinction['V'],r_v)
mw_extinction['UVW1'] = calculated_extinctions['1']
mw_extinction['UVW2'] = calculated_extinctions['2']
mw_extinction['UVM2'] = calculated_extinctions['M']

# Remove extinction correction if making data table for paper:
#for key in mw_extinction: 
#    mw_extinction[key] = 0

# Host extinction:

host_extinction = {}
extinction_translator = {'Clear':'R', 'o': 'R', 'UVW1': '1', 'UVW2': '2', 'UVM2': 'M'}
#ebv = 0
#ebv = 0.051 # From Na lines in low-res spectrum, Stefano reports E(B-V)<0.051
#ebv = 0.05 # From Na lines in high-res spectrum, Tartaglia+ report E(B-V)=0.09 \pm 0.04
av_host = ebv * r_v
  
calculated_extinctions = calc_extinction(av_host,r_v)

for key in mw_extinction:
    if key in extinction_translator:
        filt = extinction_translator[key]
    else:
        filt = key
    host_extinction[key] = calculated_extinctions[filt]
    
# Total Extinction:
    
extinction = {}
for key in mw_extinction:
    extinction[key] = mw_extinction[key]+host_extinction[key]

# Data:

lcogt1m_data = get_lcogt_data()
lcogt1m_data = correct_for_extinction(lcogt1m_data,extinction)
export_data(lcogt1m_data,'photometry_lcogt_1m.txt')
#lcogt1m_data = bin_daily(lcogt1m_data)

atel_data = get_atel_data()
# ATLAS new reductions:
atel_data[atel_data['source']=='ATLAS']['mag'][0] = 15.54
atel_data[atel_data['source']=='ATLAS']['mage'][0] = 0.14
atel_data[atel_data['source']=='ATLAS']['mag'][1] = 15.45
atel_data[atel_data['source']=='ATLAS']['mage'][1] = 0.08
atel_data = correct_for_extinction(atel_data,extinction)
#atel_data = bin_daily(atel_data)

swift_data = get_swift_data_pb()
swift_data = correct_for_extinction(swift_data,extinction)
export_data(swift_data,'photometry_swift.txt')
#swift_data = bin_daily(swift_data)

sv_data = get_sv_data()

w = np.array(['EFOSC' in r for r in sv_data['instr']]) 
pessto_data = sv_data[w]
#pessto_data = bin_daily(pessto_data)

w = np.array(['PROMPT' in r for r in sv_data['instr']]) 
dlt40_data = sv_data[w]
dlt40_data.remove_column('filter')
dlt40_data.add_column(Column(np.tile('Clear',len(dlt40_data)),name='filter'))
#dlt40_data = bin_daily(dlt40_data)

# For amateur data, add conservative 0.2 color-term error - actually not:
monard_data = {'jd':[2457653.393218], 'filter':['Clear'], 'mag':[15.354], 'mage':[np.sqrt(0.0613**2+0.0**2)]}
monard_data = correct_for_extinction(monard_data,extinction)

# For amateur data, add conservative 0.2 color-term error - actually not:
nicholls_data = {'jd':[2457652.906250,2457653.916667], 'filter':['Clear','Clear'], 'mag':[14.7939,15.4629], 'mage':[0.073630383,0.076565803]}
nicholls_data = correct_for_extinction(nicholls_data,extinction)

sanchez_unbinned_data = {'jd':[2457652.8076435188], 'filter':['Clear'], 'mag':[14.7109], 'mage':[0.353254221]}
sanchez_binned_data = {'jd':[2457652.8158410494], 'filter':['Clear'], 'mag':[14.8129], 'mage':[0.313319308]}
sanchez_data = {}
sanchez_data['jd'] = [0.5 * (sanchez_unbinned_data['jd'][0]+sanchez_binned_data['jd'][0])]
sanchez_data['mag'] = [0.5 * (sanchez_unbinned_data['mag'][0]+sanchez_binned_data['mag'][0])]
sanchez_data['mage'] = [0.5 * np.sqrt((sanchez_unbinned_data['mage'][0]**2+sanchez_binned_data['mage'][0]**2))]
sanchez_data['filter'] = ['Clear']
sanchez_data = correct_for_extinction(sanchez_data,extinction)

plot_lc()
stop

# Prepare data for Matlab fitting (Nakar & Piro 2014):

matlabdata = {}
w = lcogt1m_data['filter'] == 'V'
for key in lcogt1m_data:
    matlabdata[key] = lcogt1m_data[key][w]
w = np.array([s in ['Buso&Otero','ATLAS'] for s in atel_data['source']])
for key in matlabdata:
    matlabdata[key] = np.append(matlabdata[key],atel_data[key][w])
#for key in matlabdata:
#    matlabdata[key] = np.append(matlabdata[key],sanchez_data[key])
#for key in matlabdata:
#    matlabdata[key] = np.append(matlabdata[key],nicholls_data[key][0])
#for key in matlabdata:
#    matlabdata[key] = np.append(matlabdata[key],monard_data[key])
w = atel_data['filter'] == 'V'
for key in matlabdata:
    matlabdata[key] = np.append(matlabdata[key],atel_data[key][w])
Table.write(Table(matlabdata),'SN2016gkg_formatlab_e%s.txt' %ebv,format='ascii')

#########

# Model fitting:

fitdata = lcogt1m_data

#w = np.array([f in ['B','V','g','r','i','o'] for f in atel_data['filter']])
w = np.array([s in ['Buso&Otero','ASAS-SN','ATLAS','LCOGT'] for s in atel_data['source']])
for key in fitdata:
    fitdata[key] = np.append(fitdata[key],atel_data[key][w])

for key in fitdata:
    if key=='filter':
        fitdata[key] = np.append(fitdata[key],['swift,'+f for f in swift_data[key]])
    else:
        fitdata[key] = np.append(fitdata[key],swift_data[key])
for key in fitdata:
    fitdata[key] = np.append(fitdata[key],sanchez_data[key])
for key in fitdata:
    fitdata[key] = np.append(fitdata[key],nicholls_data[key][0])
for key in fitdata:
    fitdata[key] = np.append(fitdata[key],monard_data[key])
#for key in fitdata:
#    fitdata[key] = np.append(fitdata[key],pessto_data[key])

if velocity == 'fixed':
    ndim = 3
else:
    ndim = 4  

nwalkers = ndim*4

nw = datetime.datetime.now()
nw_str = str(nw).replace('-','').replace(' ','_').replace(':','').split('.')[0]
if fit_ni:
    postfix = '_ni'
else:
    postfix = ''
summary_filename = 'mcmcresults_' + model + '_' + velocity + '_' + weighting + '_' + str(n) + 'e' + str(ebv) + '_' + postfix + '_' + nw_str + '.txt'

# Arrange all the data and get the bandpasses:

swift_bp_path = '/Users/iair/data/UVOTfilters_Breeveld2011/'
excluded_filters = []

distinctfilters = list(set(fitdata['filter']))

all_t = {}
all_m = {}
all_merr = {}
bandpasses = {}
filternames = []

for filtname in distinctfilters:
    if filtname not in excluded_filters:
        if filtname in ['g','r','i']: 
            newfiltname = 'sdss,'+filtname
        else:
            newfiltname = filtname
        filternames.append(newfiltname)
        if 'swift' in newfiltname:
            bandpasses[newfiltname] = pysynphot.FileBandpass(swift_bp_path + newfiltname.split(',')[1]+'_UVOT.txt')
        elif newfiltname == 'o':
            bandpasses[newfiltname] = pysynphot.FileBandpass('orange_in_aa.dat')
        elif newfiltname == 'Clear':
            bandpasses[newfiltname] = pysynphot.ObsBandpass('sdss,r')
        else:
            bandpasses[newfiltname] = pysynphot.ObsBandpass(newfiltname)
        w =((fitdata['filter'] == filtname) & ((fitdata['jd']-t0)/(1+redshift)<tend_fit))
        all_t[newfiltname] = (fitdata['jd'][w]-t0)/(1+redshift)
        all_m[newfiltname] = fitdata['mag'][w]
        all_merr[newfiltname] = fitdata['mage'][w]

loadfromfile = ''
#loadfromfile = ['20161019_124541','fixedv_20161019_092700']
#loadfromfile = ['fixedv_20161019_092700']
#loadfromfile = ['sapir_free_normal_3_20161025_001858']
#loadfromfile = ['piro_fixed_early_3_20161102_233119']
#loadfromfile = ['sapir_free_early_1.5_20161108_224834','sapir_free_early_3_20161108_153529','piro_fixed_early_1.5_20161109_081425']
#loadfromfile = ['piro_free_early_3_20161110_142150','sapir_free_early_1.5_20161111_075141','sapir_free_early_3_20161110_102624'] 
#loadfromfile = ['piro_free_early_1.5e0.09__20161211_000504','sapir_free_early_1.5e0.09__20161211_211620','sapir_free_early_3e0.09__20161212_085116'] 

if loadfromfile:
    
    results_line = []
    for fn in loadfromfile:
        
        d = Table.read('mcmcresults_'+fn+'.txt',format='ascii')
        results_line.append(list(d[0].as_void()))
        
        if len(d.colnames) == 9:
            cornerdata = Table.read('mcmcchain_'+fn+'.txt',format='ascii',names=['R_e','M_e','t_offset'])    
        else:
            cornerdata = Table.read('mcmcchain_'+fn+'.txt',format='ascii',names=['R_e','M_e','v_e','t_offset'])

else:

    # initialize walkers [R_e,M_e,v_e,t_offset]
    if velocity == 'fixed':
        initial_guess = [0.4,0.02,-0.1]
        initial_spread = [0.2,0.002,0.2]
    else:
        initial_guess = [0.4,0.02,1.7,-0.1]
        initial_spread = [0.2,0.002,0.2,0.2]
    starting_positions = np.random.randn(nwalkers, ndim)*initial_spread+initial_guess
    # for t_offset, the starting position should be <=0
    for wlkr in range(nwalkers):
        starting_positions[wlkr][len(starting_positions[wlkr])-1] = np.min([starting_positions[wlkr][len(starting_positions[wlkr])-1],0])
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[filternames, bandpasses, all_t, all_m, all_merr, distance])
    pos, prob, state = sampler.run_mcmc(starting_positions, 200)
    #fig, ax = plt.subplots(4, sharex=True)
    #for i in range(4):
    #    ax[i].plot(sampler.chain[:, :, i].T, '-k', alpha=0.2)
    
    sampler.reset()
    pos, prob, state = sampler.run_mcmc(pos, nsteps)
    #print datetime.datetime.now()-n
    #fig, ax = plt.subplots(2, sharex=True)
    #for i in range(2):
    #    ax[i].plot(sampler.chain[:, :, i].T, '-k', alpha=0.2)
    
    chain_filename = 'mcmcchain_' + model + '_' + velocity + '_' + weighting + '_' + str(n) + 'e' + str(ebv) + '_' + postfix + '_' + nw_str + '.txt'
    if velocity == 'fixed':
        chain_table = Table(sampler.flatchain, names=['R_e', 'M_e', 't_offset'])
    else:
        chain_table = Table(sampler.flatchain, names=['R_e', 'M_e', 'v_e', 't_offset'])
    chain_table.write(chain_filename,format='ascii.fast_no_header')
    
    #corner.corner(sampler.flatchain, labels=['R_e', 'v_e', 'm_e', 't_offset'])
    #triangle_plot.plot(sampler.flatchain.T,chain_filename.replace('.txt','.pdf'),labels=['Temperature', 'Radius'], scientific_notation=[False,True])
    
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    
    results = zip(*np.percentile(samples, [16, 50, 84],axis=0))
    if velocity == 'fixed':
        results_line = np.reshape(results,[1,9])
    else:
        results_line = np.reshape(results,[1,12])
    print(results_line)
    print('')
    
    if velocity == 'fixed':    
        s = '{} {} {} {} {} {} {} {} {}\n'.format(*results_line[0])
    else:
        s = '{} {} {} {} {} {} {} {} {} {} {} {}\n'.format(*results_line[0])
    with open(summary_filename,'w') as f: 
        f.write(s)
        
    results_line = [results_line[0]]
    cornerdata = chain_table
    
    #try:
    #    fig = corner.corner(samples, plot_contours=False, labels=['$R_e$ [$10^{13}$ cm]', '$v_e$ [$10^4$ km/s]', '$m_e$ [$M_{\odot}$]', '$t_{offset}$ [days]'])
    #    fig.savefig(figpath+chain_filename.replace('.txt','.pdf'))
    #except Exception as e:
    #    print e

best_re = []
best_me = []
best_ve = []
best_toffset = []
for rn in results_line:
    best_re.append(rn[0:3])
    best_me.append(rn[3:6])
    if len(rn) == 9:
        best_ve.append([1.5,1.5,1.5])
        best_toffset.append(rn[6:9])
    else:
        best_ve.append(rn[6:9])
        best_toffset.append(rn[9:12])
fitparams = [[best_re[i][1]*1e13 for i in range(len(results_line))],[best_me[i][1] for i in range(len(results_line))],[best_ve[i][1]*1e9 for i in range(len(results_line))],[best_toffset[i][1] for i in range(len(results_line))]]
plot_lc(*fitparams,bersten=False)


cornerdata['R_e'] = cornerdata['R_e']*10.0
cornerdata['M_e'] = cornerdata['M_e']*100.0
if velocity == 'fixed':
    labels = ['$R_e$ [$10^{12}\\,$cm]', '$M_e$ [$10^{-2}\\,M_{\odot}$]', '$t_{\\rm offset}$ [days]']
else:
    labels = ['$R_e$ [$10^{12}\\,$cm]', '$M_e$ [$10^{-2}\\,M_{\odot}$]', '$v_e$ [$10^4\\,$km$\\,s^{-1}$]', '$t_{\\rm offset}$ [days]']
fig,axs = beautifulcorner.corner(cornerdata,labels=labels)
if len(loadfromfile)>=1:
    fig.savefig(figpath+'corner_' + loadfromfile[0] + '.pdf', bbox_inches='tight')
else:
    fig.savefig(figpath+'corner_' + summary_filename.replace('mcmcresults','corner').replace('.txt','') + '.pdf', bbox_inches='tight')

#rsun = astropy.constants.R_sun.to('cm').value
#results = {}
#results['params'] = ['$R_e$','','$M_e$','$v_e$','$t_{\\rm offset}$']
#results['values'] = ['$%0.2f{\\times}10^{12}$' %(best_re[1]*10), '$%0.1f$' %(best_re[1]*1e13/rsun), '$%0.2f{\\times}10^{-2}$' %(best_me[1]*100), '$%0.2f{\\times}10^{4}$' %(best_ve[1]),'$%0.3f$' %(best_toffset[1])]
#results['bounds'] = ['$\left(%0.2f-%0.2f\\right){\\times}10^{12}$' %(best_re[0]*10,best_re[2]*10), '$\left(%0.1f-%0.1f\\right)$' %(best_re[0]*1e13/rsun,best_re[2]*1e13/rsun), '$\left(%0.2f-%0.2f\\right){\\times}10^{-2}$' %(best_me[0]*100,best_me[2]*100), '$\left(%0.2f-%0.2f\\right){\\times}10^{4}$' %(best_ve[0],best_ve[2]),'$\left(%0.3f\\right)-\left(%0.3f\\right)$' %(best_toffset[0],best_toffset[2])]
#results['units'] = ['cm','$\\Rsun$','$\Msun$','km\,s$^{-1}$','days']

## Matlab Nakar & Piro (2014) results [ASSUMES v_e=1e9]:
#
## Lpeak, Tpeak, Lmin, m_e [10^-3 Msun], r_e [10^13 cm], r_core [10^11 cm]
## {SN2016gkg} & {$1.69^{+0.06}_{-0.06}$} & {$1.33^{+0.02}_{-0.02}$} & {$0.38^{+0.00}_{-0.38}$} & {$8.84^{+0.27}_{-0.26}$} & {$2.27\pm0.03$} & {$<167.55$}\tabularnewline
#
#np14M_e = [8.84e-3, -0.26e-3, 0.27e-3]
#np14R_e = [2.27e13, -0.03e13, 0.03e13]
#best_v9 = [best_ve[1]/1.0e9, (best_ve[1]-best_ve[0])/1.0e9]
#corrected_np14M_e = [np14M_e[0]*(best_v9[0])]
#corrected_np14R_e = [np14R_e[0]/((best_v9[0])**(2))]
#
## Errors:
#corrected_np14M_e.append(np.sqrt((np14M_e[2]*best_v9[0])**2+(np14M_e[0]*best_v9[1])**2))
#corrected_np14R_e.append(np.sqrt((np14R_e[2]/(best_v9[0]**2))**2+(2*np14R_e[0]*best_v9[1]/(best_v9[0]**3))**2))
#
#results_np14 = {}
#results_np14['params'] = ['$R_e$','','$M_e$']
#results_np14['values'] = ['$%0.2f{\\times}10^{12}v_{e,4}^{-2}$' %(np14R_e[0]/1e12), '$%0.1fv_{e,4}^{-2}$' %(np14R_e[0]/rsun), '$%0.2f{\\times}10^{-2}v_{e,4}$' %(np14M_e[0]*100)]
#results_np14['bounds'] = ['$\left(%0.2f-%0.2f\\right){\\times}10^{12}v_{e,4}^{-2}$' %(np14R_e[0]/1e12+np14R_e[1]/1e12,np14R_e[0]/1e12+np14R_e[2]/1e12), '$\left(%0.1f-%0.1f\\right)v_{e,4}^{-2}$' %((np14R_e[0]+np14R_e[1])/rsun,(np14R_e[0]+np14R_e[2])/rsun), '$\left(%0.2f-%0.2f\\right){\\times}10^{-2}v_{e,4}$' %((np14M_e[0]+np14M_e[1])*100,(np14M_e[0]+np14M_e[2])*100)]
#results_np14['units'] = ['cm','$\\Rsun$','$\Msun$']
#
#resultstable = Table([results[field]+results_np14[field] for field in ['params','values','bounds','units']], names=['Parameter','Value','$68\%$ CB','Unit'])
#resultstable.write('results.tex',format='latex')
#
#print 'Explosion time:' 
#print [t for t in t0+np.array(best_toffset)]

