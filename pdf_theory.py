from scipy import *
import numpy as np
from emcee.utils import MPIPool 
import sys, itertools
import emcee
import os
from astropy.io import fits
import scipy.ndimage as snd

#zidx=int(sys.argv([1])) ### index for redshift, go from 0 to 5

smooth = lambda kmap, sigma: snd.filters.gaussian_filter(kmap.astype(float),sigma,mode='constant')

map_dir = '/tigress/jialiu/will/convergence_6redshifts/'
out_dir = '/tigress/jialiu/PDFs/cora/'

mapgen = lambda z, r, mnu: fits.open(map_dir+'convergence_6redshifts_mnv0.%i0000_om0.30000_As2.1000/Maps%02d/WLconv_z%.2f_%04dr.fits'%(mnu, z*10, z, r))[0].data ## mnu=1 or 0

thetaG_arr = array([5, 10, 15]) * 512./(3.5*60) #arcmin * pix_per_arcmin

z_arr = concatenate([arange(0.5, 3, 0.5),[1100.,]])
#sigmakappa_arr = array([ [std(smooth(mapgen(iz, 1, 0), thetaG)) 
                  #for thetaG in thetaG_arr] for iz in z_arr])

sigmakappa_arr = array([[ 0.00320448,  0.0022913 ,  0.00168369],
                        [ 0.0051601 ,  0.00390173,  0.00299778],
                        [ 0.00667598,  0.00507506,  0.0039651 ],
                        [ 0.00764257,  0.00567065,  0.00433421],
                        [ 0.00856384,  0.00624893,  0.0046596 ],
                        [ 0.01968398,  0.01379884,  0.00995027]])

##### pre-defined bin edges, 100 linear bins between -5, 5 sigma_kappa
binedges_fun = lambda sigmak: np.linspace(-5*sigmak, 5*sigmak, 101)
binedges = array([map(binedges_fun,isigmakappa_arr) for isigmakappa_arr in sigmakappa_arr])

#### smooth function, for all smoothing bins
def smooth_map (r, mnu=0, zidx=0):
    imap = mapgen(z_arr[zidx], r, mnu)
    imap_smooth = array([smooth(imap, thetaG) for thetaG in thetaG_arr])
    hist_arr = array([histogram(imap_smooth[i], bins=binedges[zidx, i])[0] 
                      for i in range(len(thetaG_arr))])
    return hist_arr
    
pool=MPIPool()
if not pool.is_master():
    pool.wait()
    sys.exit(0)

#for imnu in (0,1):
    #for izidx in range(6):
imnu, izidx = 0, 0
print 'Mnu, z:', imnu, z_arr[izidx]
out = array(pool.map(lambda p: smooth_map(p, mnu=imnu, zidx=izidx), range(1,1001)))
for j in range(len(thetaG_arr)):
    save(out_dir+'PDFs_Mnu0.%i_z%.1f_smooth%02d.npy'%(imnu, z_arr[izidx], thetaG_arr[j]), out[:,j,:])

pool.close()

#out_massless = array(pool.map(lambda p: smooth_map(p, mnu=1), range(1,10000)))
### output shape=(9999, 3, 100) for 9999 realizations, 3 smoothing scale, 100 bins
