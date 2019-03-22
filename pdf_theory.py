from scipy import *
import numpy as np
from emcee.utils import MPIPool 
import sys, itertools
import emcee
import os
from astropy.io import fits
import scipy.ndimage as snd

smooth = lambda kmap, sigma: snd.filters.gaussian_filter(kmap.astype(float),sigma,mode='constant')

map_dir = '/tigress/jialiu/will/convergence_6redshifts/'

mapgen = lambda z, r, mnu: fits.open(map_dir+'convergence_6redshifts_mnv0.%i0000_om0.30000_As2.1000/Maps%02d/WLconv_z%.2f_%04dr.fits'%(mnu, z*10, z, r))[0].data ## mnu=1 or 0

thetaG_arr = array([10,15,20]) * 512./(3.5*60) #arcmin * pix_per_arcmin

z_arr = concatenate([arange(0.5, 3, 0.5),[1100.,]])
#sigmakappa_arr = array([ [std(smooth(mapgen(iz, 1, 0), thetaG)) 
                  for thetaG in thetaG_arr] for iz in z_arr])

sigmakappa_arr = array([[ 0.00320448,  0.0022913 ,  0.00168369],
                        [ 0.0051601 ,  0.00390173,  0.00299778],
                        [ 0.00667598,  0.00507506,  0.0039651 ],
                        [ 0.00764257,  0.00567065,  0.00433421],
                        [ 0.00856384,  0.00624893,  0.0046596 ],
                        [ 0.01968398,  0.01379884,  0.00995027]])
binedges = lambda 
def smooth_map (mnu, z, r, thetaG):
    imap = mapgen(z, r, mnu)
    imap_smooth = WLanalysis.smooth(imap, thetaG)
    

if not plot_only:
    pool=MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
