from scipy import *
import numpy as np
import WLanalysis
from emcee.utils import MPIPool 
import sys, itertools
import emcee
import os

local=0
tightball = 0
add_2dpdf = 0
plot_only = 0
single_z = 0
test_cross = 0
very_wide = 0
upload_dropbox = 1
upload_MCMC=1
bloc_cov = 0

Nk='10k' # '5ka', '5kb'
Nmin=1000 #500###### minimum counts in that bin to get included in PDF calculation
Nmin2=20
Nchain = 1000
iscale = 1 ## rescale the PDF so it has similar magnitude as the power spectrum

#Nmin_scale_arr = [[iNmin, iscale] for iscale in (1,1e-12, 1e-14) 
                #for iNmin in (500, 1000, 5000, 100) ]

try:
    Nk = str(sys.argv[1])
    single_z=int(sys.argv[2])
    #Nmin,iscale=Nmin_scale_arr [int(sys.argv[2])]
    
except Exception:
    pass

collapse=''#'collapsed'#
np.random.seed(10026)#

testfn = collapse+'Mar6_%s_ell600_%s_Nchain%i_%s_Nmin%i'%(['tomo','z1'][single_z],['wideP0','tightball'][tightball],Nchain,Nk,Nmin)#''#

if very_wide:
    testfn = collapse+'Sep5_%s_fullcov_%s_Nchain%i_%s'%(['tomo','z1'][single_z],['verywideP0','tightball'][tightball],Nchain,Nk)#''#
#testfn = collapse+'Aug16_R_Nmin%s_Nmin2%s_Nchain%i_%s'%(Nmin,Nmin2,Nchain,Nk)#''#
Nmin*=iscale

z_arr = arange(0.5,3,0.5)

#####################################
######## set up folders #############
#####################################

######## stampede2
like_dir='/scratch/02977/jialiu/peakaboo/likelihood_tomo/'
stats_dir = '/scratch/02977/jialiu/peakaboo/stats_tomo/'
peakaboo_dir='/scratch/02977/jialiu/peakaboo/'

if single_z:
    z_arr = [1.0,]
    like_dir='/scratch/02977/jialiu/peakaboo/likelihood_z1/'
    stats_dir = '/scratch/02977/jialiu/peakaboo/stats_z1/'

#ebcov_dir = stats_dir+'Om0.29997_As2.10000_mva0.00000_mvb0.00000_mvc0.00000_h0.70000_Ode0.69995/1024b512/box5/output_eb_5000_s4/seed0/'

ebcov_dir = stats_dir+'box5/'

Nz = len(z_arr)

######### local
if local:
    #stats_dir = '/Users/jia/Dropbox/weaklensing/PDF/'
    #ebcov_dir = stats_dir+'box5/output_eb_5000_s4/seed0/'
    ebcov_dir = '/Users/jia/Dropbox/weaklensing/PDF/box5/'
    stats_dir = '/Users/jia/Dropbox/weaklensing/PDF/'
    peakaboo_dir = '/Users/jia/Dropbox/weaklensing/PDF/'
    like_dir = peakaboo_dir+'likelihood'
    upload_dropbox = 0
    
eb_dir = stats_dir+'stats_avg/output_eb_5000_s4/'
#eb1k_dir = stats_dir+'stats_avg_1k/output_eb_5000_s4/'
#

params = genfromtxt(peakaboo_dir+'cosmo_params_all.txt',usecols=[2,3,4])
#####################################
##### initiate avg statistics #######
#####################################

###### PS shape:(15, 101, 20)

# auto's only
psIauto = array( [load(eb_dir+'ALL_igalXigal_z{0}_z{0}_{1}.npy'.format(iz,Nk)) for iz in z_arr])

##### 1d PDF shape:(5, 101, 27)
pdf1dN = iscale*array( [load(eb_dir+'ALL_gal_pdf_z{0}_sg1.0_{1}.npy'.format(iz,Nk)) for iz in z_arr])

#### 2d PDF shape:(10, 101, 27, 27)
if add_2dpdf:
    pdf2dN = array( [load(eb_dir+'ALL_galXgal_2dpdf_z{0}_z{1}_sg1.0_{2}.npy'.format(z_arr[i],z_arr[j],Nk)) 
                    for i in range(Nz) for j in range(i+1,Nz)])
    ##pdf2dN1ks = array( [[load(eb1k_dir+'ALL_galXgal_2dpdf_z{0}_z{1}_sg1.0_1k{2}.npy'.format(z_arr[i],z_arr[j], ik)) 
                    ##for i in range(Nz) for j in range(i+1,Nz)] for ik in range(10)])

    ########### test collapsed 1d PDF from 2d shape:(5, 101, 27), mean
    if collapse=='collapsed':
        pdf1dNb = array([sum(pdf2dN[i],axis=-1) for i in [0,4,7,9] ] + [sum(pdf2dN[-1],axis=-2)])

#####################################
###### covariances stats ############
#####################################

covIgen = lambda ips_cov:mat(cov(ips_cov,rowvar=0)*12.25/2e4).I
    
##### PS shape:(101,100)
psIauto_flat = swapaxes(psIauto,0,1).reshape(101,-1) 
psNauto_cov = swapaxes(array( [load(ebcov_dir+'ALL_galXgal_z{0}_z{0}.npy'.format(z_arr[i]))
                           for i in range(Nz)]),0,1).reshape(10000,-1)
covIpsNauto = covIgen(psNauto_cov)

###### PDF 1D
idxt=where(pdf1dN[:,1]>Nmin)

pdf1dN_flat= swapaxes(pdf1dN[idxt[0],:,idxt[1]],0,1).reshape(101,-1) 
##pdf1dN1k_flat = array([swapaxes(ips[idxt[0],:,idxt[1]],0,1).reshape(101,-1) for ips in pdf1dN1ks])

pdf1dN_cov = iscale*swapaxes(array( [load(ebcov_dir+'ALL_gal_pdf_z{0}_sg1.0.npy'.format(iz)) for iz in z_arr])[idxt[0],:,idxt[1]],0,1).reshape(10000,-1)
covIpdf1dN = covIgen(pdf1dN_cov)

###### combined ps + pdf, for both auto and cross

def covIgen_diag (psN_cov, N=pdf1dN_flat.shape[-1]): ##### test zero out off-diag terms
    covfull = cov(psN_cov,rowvar=0)*12.25/2e4
    ######### last ? are pdf
    covdiag = zeros(shape=covfull.shape)
    covdiag[:-N,:-N]=1
    covdiag[-N:,-N:]=1
    covdiag*=covfull
    return mat(covdiag).I


comb_auto_flat = concatenate([psIauto_flat, pdf1dN_flat], axis=-1)
#comb_cros_flat = concatenate([psI_flat, pdf1dN_flat], axis=-1)

comb_cov_auto = concatenate([psNauto_cov,pdf1dN_cov],axis=-1)
#comb_cov_cros = concatenate([psN_cov,pdf1dN_cov],axis=-1)

covIcomb_auto = covIgen(comb_cov_auto)
#covIcomb_cros = covIgen(comb_cov_cros)

#covIcomb_auto = covIgen_diag(comb_cov_auto)
#covIcomb_cros = covIgen_diag(comb_cov_cros)

###### PDF 2D
if add_2dpdf:
    idxt2=where(pdf2dN[:,5]>Nmin2)

    pdf2dN_flat= swapaxes(pdf2dN,0,1)[:,idxt2[0],idxt2[1],idxt2[2]]
    #pdf2dN1k_flat= array([swapaxes(ips,0,1)[:,idxt2[0],idxt2[1],idxt2[2]] for ips in pdf2dN1ks])

    pdf2dN_cov = swapaxes(array( [load(ebcov_dir+'ALL_galXgal_2dpdf_z{0}_z{1}_sg1.0.npy'.format(z_arr[i],z_arr[j]))
                                for i in range(Nz) for j in range(i+1,Nz)]),0,1)[:,idxt2[0],idxt2[1],idxt2[2]].reshape(10000,-1)

    covpdf2dN = cov(pdf2dN_cov,rowvar=0)*12.25/2e4
    covIpdf2dN = mat(covpdf2dN).I

    ############### test collapsed 1d PDF from 2d, covariance
    if collapse=='collapsed':
        pdf1dN_cov = array([sum(load(ebcov_dir+'ALL_galXgal_2dpdf_z{0}_z{1}_sg1.0.npy'.format(z_arr[i],z_arr[i+1])),axis=-1) for i in range(4) ] + [sum(load(ebcov_dir+'ALL_galXgal_2dpdf_z2.0_z2.5_sg1.0.npy'),axis=-2)])[idxt[0],:,idxt[1]].T
        covpdf1dN = cov(pdf1dN_cov,rowvar=0)*12.25/2e4
        covIpdf1dN = mat(covpdf1dN).I


#####################################
###### build emulator ###############
#####################################

fidu_params = array([0.1,0.3,2.1])

#obss = [psI_flat[1], pdf1dN_flat[1], pdf2dN_flat[1]]
#covIs = [covIpsN, covIpdf1dN, covIpdf2dN]

#obss = [psIauto_flat[1], psI_flat[1], pdf1dN_flat[1], comb_auto_flat[1],comb_cros_flat[1]]
#covIs = [covIpsNauto, covIpsN, covIpdf1dN, covIcomb_auto, covIcomb_cros]

#emusingle = [WLanalysis.buildInterpolator(array(istats)[1:], params[1:], function='GP') 
             #for istats in [psIauto_flat, psI_flat, pdf1dN_flat]] #comb_auto_flat,comb_cros_flat]]
#emucomb_auto = lambda p: concatenate([emusingle[0](p),emusingle[2](p)])
#emucomb_cross = lambda p: concatenate([emusingle[1](p),emusingle[2](p)])
#emulators= emusingle+[emucomb_auto,emucomb_cross]

#emulators = [WLanalysis.buildInterpolator(array(istats)[1:], params[1:], function='GP') 
#             for istats in [psIauto_flat, psI_flat, pdf1dN_flat, comb_auto_flat,comb_cros_flat]]

if test_cross:
    ######### cross
    psIcross = array( [load(eb_dir+'ALL_igalXigal_z{0}_z{1}_{2}.npy'.format(z_arr[i],z_arr[j],Nk))
              for i in range(Nz-1) for j in range(i+1,Nz)])
    psIcross_flat = swapaxes(psIcross,0,1).reshape(101,-1) 
    psNcross_cov = swapaxes(array( [load(ebcov_dir+'ALL_galXgal_z{0}_z{1}.npy'.format(z_arr[i],z_arr[j]))
                           for i in range(Nz-1) for j in range(i+1,Nz)]),0,1).reshape(10000,-1)
    
    covIpsNcross = covIgen(psNcross_cov)
    
    ######### auto + cross 
    #psI = array( [load(eb_dir+'ALL_igalXigal_z{0}_z{1}_{2}.npy'.format(z_arr[i],z_arr[j],Nk))
              #for i in range(Nz) for j in range(i,Nz)])
    #psI_flat = swapaxes(psI,0,1).reshape(101,-1) 
    #psN_cov = swapaxes(array( [load(ebcov_dir+'ALL_galXgal_z{0}_z{1}.npy'.format(z_arr[i],z_arr[j]))
                          #for i in range(Nz) for j in range(i,Nz)]),0,1).reshape(10000,-1)
    #covIpsN = covIgen(psN_cov)
    
    psI_flat = concatenate([psIauto_flat, psIcross_flat], axis=-1)
    psN_cov = concatenate([psNauto_cov,psNcross_cov],axis=-1)
    covIpsN = covIgen(psN_cov)
    
    obss = [psIauto_flat[1],psIcross_flat[1],psI_flat[1]]
    covIs = [covIpsNauto, covIpsNcross, covIpsN]
    emulators = [WLanalysis.buildInterpolator(array(istats)[1:], params[1:], function='GP') 
                for istats in [psIauto_flat,  psIcross_flat, psI_flat]]
    
else:
    obss = [psIauto_flat[1], pdf1dN_flat[1], comb_auto_flat[1]]
    covIs = [covIpsNauto, covIpdf1dN, covIcomb_auto]

    emusingle = [WLanalysis.buildInterpolator(array(istats)[1:], params[1:], function='GP') 
                for istats in [psIauto_flat,  pdf1dN_flat]] #comb_auto_flat,comb_cros_flat]]
    emucomb_auto = lambda p: concatenate([emusingle[0](p),emusingle[1](p)])
    emulators= emusingle+[emucomb_auto,]

if add_2dpdf:
    obss += [pdf2dN_flat,]
    covIs +=[covIpdf2dN,]
    emulators += [WLanalysis.buildInterpolator(array(pdf2dN_flat)[1:], params[1:], function='GP'),]

#emulators = [WLanalysis.buildInterpolator(array(istats)[1:], params[1:], function='GP') 
            #for istats in [psIauto_flat, pdf1dN_flat, comb_auto_flat]]

#########
rDH = [ float((1e4-len(covI)-2.0)/9999.0) for covI in covIs] ## 

############# likelihood #######
from scipy.misc import factorial


def lnprob_gaussian(p,jjj):
    '''log likelihood '''
    if p[0]<0: ####### force neutrino mass to be positive
        return -np.inf
    diff = emulators[jjj](p)-obss[jjj]
    return float(-0.5*mat(diff)*covIs[jjj]*mat(diff).T)*rDH[jjj]

def lnprob_block(p,jjj=0):
    '''block likelihood 
    '''
    return lnprob_gaussian(p,0)+lnprob_gaussian(p,1)

def lnprob_poisson(p,jjj=1):
    if p[0]<0: ####### force neutrino mass to be positive
        return -np.inf
    mu = emulators[jjj](p)
    n = obss[jjj]
    return sum(n*log(mu)-mu-log(factorial(n)))


lnprob = lnprob_gaussian

#fn_arr = ['psAuto','psCross','pdf1d','combAuto','combCross']
if not test_cross:
    fn_arr = ['psAuto','pdf1d','combAuto']#,'pdf2d','comb2d']
else:
    fn_arr = ['psAuto','psCross','Auto+cross']
if add_2dpdf:
    fn_arr += ['pdf2d',]
if not plot_only:
    pool=MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    print Nk, int(Nmin/iscale), iscale, testfn

    nwalkers=544
    ndim=3
    #p0 = (array([ (rand(nwalkers, ndim) -0.5) * array([1, 0.3, 0.3]) + 1]) * fidu_params).reshape(-1,3)
    
    ########## wide
    p0_ranges=array([[0,0.4],[0.25,0.35],[1.6,2.6]])
    ########### very wide
    if very_wide:
        p0_ranges=array([[0.0,1.0],[0.1,0.6],[1.2,3.0]])
    p0=rand(nwalkers,ndim)*(p0_ranges[:,1]-p0_ranges[:,0]).reshape(1,3)+p0_ranges[:,0].reshape(1,3)
    ########### tight ball
    if tightball:
        p0 = (array([ (rand(nwalkers, ndim) -0.5) * 1e-2 * array([1, 0.3, 0.3]) + 1]) * fidu_params).reshape(-1,3)
    
    #print 'rDH',rDH
    #for i in range(len(covIs)):
        #print fn_arr[i]
        #print 'cov shape',covIs[i].shape
        #print 'stats shape',[psIauto_flat,  pdf1dN_flat, comb_auto_flat][i].shape
        
    i=0
    print fn_arr[i], obss[i].shape
    ifn = like_dir+'MC_%s_%s.npy'%(fn_arr[i],testfn)
    if not os.path.isfile(ifn):
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[i,], pool=pool)
        pos, prob, state = sampler.run_mcmc(p0, 100)
        sampler.reset()
        sampler.run_mcmc(pos, Nchain)
        save(ifn, sampler.flatchain)

    i=1
    print fn_arr[i], obss[i].shape
    ifn = like_dir+'MC_%s_%s.npy'%(fn_arr[i],testfn)
    if not os.path.isfile(ifn):
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[i,], pool=pool)
        pos, prob, state = sampler.run_mcmc(p0, 100)
        sampler.reset()
        sampler.run_mcmc(pos, Nchain)
        save(ifn, sampler.flatchain)

    i=2
    if bloc_cov:
        lnprob = lnprob_block
    print fn_arr[i], obss[i].shape
    ifn = like_dir+'MC_%s_%s.npy'%(fn_arr[i],testfn)
    if not os.path.isfile(ifn):
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[i,], pool=pool)
        pos, prob, state = sampler.run_mcmc(p0, 100)
        sampler.reset()
        sampler.run_mcmc(pos, Nchain)
        save(ifn, sampler.flatchain)

    if add_2dpdf:
        i=3
        print fn_arr[i]
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[i,], pool=pool)
        pos, prob, state = sampler.run_mcmc(p0, 100)
        sampler.reset()
        sampler.run_mcmc(pos, Nchain)
        save(like_dir+'MC_%s_%s.npy'%(fn_arr[i],testfn), sampler.flatchain)

    #i=4
    #print fn_arr[i]
    #sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[i,], pool=pool)
    #pos, prob, state = sampler.run_mcmc(p0, 100)
    #sampler.reset()
    #sampler.run_mcmc(pos, Nchain)
    #save(like_dir+'MC_%s_%s.npy'%(fn_arr[i],testfn), sampler.flatchain)

    #print 'PDF 2D'
    ##p0 = (array([ (rand(nwalkers, ndim) -0.5) * array([1, 0.3, 0.3]) + 1]) * fidu_params).reshape(-1,3)
    #sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[2,], pool=pool)
    #pos, prob, state = sampler.run_mcmc(p0, 100)
    #sampler.reset()
    #sampler.run_mcmc(pos, Nchain*10)
    #save(like_dir+'MC_%s_%s.npy'%(fn_arr[i],testfn), sampler.flatchain)

    pool.close()
    print 'done done done'

########### plotting and uploading to dropbox 
import corner
from matplotlib import pyplot
from scipy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import *

fidu_params = array([0.1,0.3,2.1])
colors=['c','y','r','b','g']
proxy=[plt.Rectangle((0,0),1,0.5,ec=icolor, fc = icolor) for icolor in colors]

stats_dir = '/scratch/02977/jialiu/peakaboo/'
#range=[[-0.1,0.45],[0.28,0.32],[1.8,2.7]]
def plotmc(chain, f=None, icolor='k',range=[[-0.1,0.5],[0.27,0.33],[1.7,2.7]]):
    chain = chain[len(chain)/3:]
    corner.corner(chain, labels=[r"$M_\nu$", r"$\Omega_m$", r"$A_s$"],
                  levels=[0.95,],color=icolor,
                   range=range,truths=fidu_params, fig=f, 
                  plot_datapoints=0, plot_density=0,
                  truth_color="k",fill_contours=0)#0.67,

#MC_ps = [load(like_dir+'MC_ps_base.npy'),]
##MC_arr = MC_ps+[load(like_dir+'MC_%s_%s.npy'%(ips,testfn)) for ips in
##               fn_arr[1:]]

#MC_arr = MC_ps + [load(like_dir+'MC_pdf1d_Aug16_tightball_R_Nmin500_Nchain500_10k.npy'),]
#MC_arr = MC_arr + [load(like_dir+'MC_%s_%s.npy'%(fn_arr[2],testfn)),]

MC_arr = [load(like_dir+'MC_%s_%s.npy'%(ips,testfn)) for ips in fn_arr]

f,ax=subplots(3,3,figsize=(6,6))
for j in range(len(MC_arr)):
    plotmc(MC_arr[j],f=f,icolor=colors[j])
ax[0,1].legend(proxy[:len(fn_arr)],fn_arr,fontsize=8)
ax[0,1].set_title(testfn)
fnfig='contour_%s.jpg'%(testfn)
fnpath=like_dir+'plots/'+fnfig
savefig(fnpath)
close()

if upload_dropbox:
    print 'uploading to dropbox'
    import os
    os.system('/work/02977/jialiu/Dropbox-Uploader/dropbox_uploader.sh upload %s %s'%(fnpath,fnfig))
    if upload_MCMC:
        for ips in fn_arr:
            ifn='MC_%s_%s.npy'%(ips,testfn)
            os.system('/work/02977/jialiu/Dropbox-Uploader/dropbox_uploader.sh upload %s %s'%(like_dir+ifn,ifn))
