from scipy import *
import numpy as np
from emcee.utils import MPIPool 
import sys, itertools
import emcee
import os
from astropy.io import fits
import scipy.ndimage as snd

#zidx=int(sys.argv([1])) ### izidx index for redshift, go from 0 to 5
imnu, izidx = int(sys.argv[1]), int(sys.argv[2])

smooth = lambda kmap, sigma: snd.filters.gaussian_filter(kmap.astype(float),sigma,mode='constant')

map_dir = '/tigress/jialiu/will/convergence_6redshifts/'
out_dir = '/tigress/jialiu/PDFs/cora/'

mapgen = lambda z, r, mnu: fits.open(map_dir+'convergence_6redshifts_mnv0.%i0000_om0.30000_As2.1000/Maps%02d/WLconv_z%.2f_%04dr.fits'%(mnu, z*10, z, r))[0].data ## mnu=1 or 0

thetaG_arcmin = array([5, 10, 15]) #arcmin * pix_per_arcmin

thetaG_arr = thetaG_arcmin * 512./(3.5*60) #arcmin * pix_per_arcmin

z_arr = concatenate([arange(0.5, 3, 0.5),[1100.,]])
sigmakappa_arr = array([ [std(smooth(mapgen(iz, 5, 0), thetaG).flatten()) 
                 for thetaG in thetaG_arr] for iz in z_arr])

sigmakappa_arr_massive = array([[ 0.00459626,  0.00308772,  0.00220381],
       [ 0.00708774,  0.00497572,  0.003762  ],
       [ 0.00927037,  0.00648126,  0.00492705],
       [ 0.01086472,  0.00743865,  0.00551859],
       [ 0.01226901,  0.00834649,  0.00608883],
       [ 0.03019188,  0.01931345,  0.01356282]])


sigmakappa_arr = array([[ 0.0047579 ,  0.00320448,  0.0022913 ],
       [ 0.00734598,  0.0051601 ,  0.00390173],
       [ 0.00954841,  0.00667598,  0.00507506],
       [ 0.01116931,  0.00764257,  0.00567065],
       [ 0.01260201,  0.00856384,  0.00624893],
       [ 0.03081554,  0.01968398,  0.01379884]])

correct_ratio = array([[0.7652519052760128, 0.8536888729600138, 0.98560815],
       [0.97673189, 1.01999364, 1.09682727],
       [1.0377963 , 1.08292562, 1.15041094],
       [1.08185203, 1.14821166, 1.24154902],
       [1.0987346 , 1.16937103, 1.27900662],
       [1.16195681, 1.24064477, 1.34414003]])

sigmakappa_arr*=correct_ratio

##### pre-defined bin edges, 100 linear bins between -5, 5 sigma_kappa
binedges_fun = lambda sigmak: np.linspace(-5*sigmak, 5*sigmak, 101)
binedges = array([map(binedges_fun,isigmakappa_arr) for isigmakappa_arr in sigmakappa_arr])

#### smooth function, for all smoothing bins
def smooth_map (r, mnu=imnu, zidx=izidx):
    imap = mapgen(z_arr[zidx], r, mnu)
    imap_smooth = array([smooth(imap, thetaG) for thetaG in thetaG_arr])
    hist_arr = array([histogram(imap_smooth[i], bins=binedges[zidx, i])[0] 
                      for i in range(len(thetaG_arr))])
    std_arr = array([std(imap) for imap in imap_smooth])
    return hist_arr, std_arr
    
pool=MPIPool()
if not pool.is_master():
    pool.wait()
    sys.exit(0)

print 'Mnu, z:', imnu, z_arr[izidx]
out_all = array(pool.map(smooth_map, range(1,10001)))
out = array([out_all[i,0] for i in range(len(out_all))])
out_std = array([out_all[i,1] for i in range(len(out_all))])
for j in range(len(thetaG_arr)):
    save(out_dir+'PDFs_Mnu0.%i_z%.1f_smooth%02d.npy'%(imnu, z_arr[izidx], thetaG_arcmin[j]), out[:,j,:])
save(out_dir+'PDFs_Mnu0.%i_z%.1f.npy'%(imnu, z_arr[izidx]), out_std)

pool.close()

#out_massless = array(pool.map(lambda p: smooth_map(p, mnu=1), range(1,10000)))
### output shape=(9999, 3, 100) for 9999 realizations, 3 smoothing scale, 100 bins

print 'DONE DONE DONE'
