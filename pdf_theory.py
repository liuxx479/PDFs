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

sigmakappa_arr_massive = array([[0.00356207, 0.00249696, 0.00195413],
       [0.0067694 , 0.00477232, 0.00374541],
       [0.00935179, 0.00661128, 0.00517522],
       [0.01143028, 0.00807851, 0.0063025 ],
       [0.01312191, 0.00925983, 0.00720139],
       [0.03459323, 0.02328923, 0.01740048]])


sigmakappa_arr = array([[0.00368035, 0.00257658, 0.00201411],
       [0.00697669, 0.00490864, 0.00384789],
       [0.00962323, 0.00678943, 0.00530901],
       [0.01175002, 0.0082883 , 0.0064598 ],
       [0.01347936, 0.0094942 , 0.00737679],
       [0.03531527, 0.02374736, 0.01773231]])

##### pre-defined bin edges, 100 linear bins between -5, 5 sigma_kappa
binedges_fun = lambda sigmak: np.linspace(-5*sigmak, 5*sigmak, 101)
binedges = array([map(binedges_fun,isigmakappa_arr) for isigmakappa_arr in sigmakappa_arr])

#### smooth function, for all smoothing bins
def smooth_map (r, mnu=imnu, zidx=izidx):
    imap = mapgen(z_arr[zidx], r, mnu)
    imap_smooth = array([smooth(imap, thetaG) for thetaG in thetaG_arr])
    hist_arr = array([histogram(imap_smooth[i], bins=binedges[zidx, i])[0] 
                      for i in range(len(thetaG_arr))])
    return hist_arr

def std_map (r, mnu=imnu, zidx=izidx):
    imap = mapgen(z_arr[zidx], r, mnu)
    imap_smooth = array([smooth(imap, thetaG) for thetaG in thetaG_arr])
    std_arr = [std(imap) for imap in imap_smooth]
    return std_arr

pool=MPIPool()
if not pool.is_master():
    pool.wait()
    sys.exit(0)

print 'Mnu, z:', imnu, z_arr[izidx]
out = array(pool.map(smooth_map, range(1,10001)))
for j in range(len(thetaG_arr)):
    print out.shape
    save(out_dir+'PDFs_Mnu0.%i_z%.1f_smooth%02d.npy'%(imnu, z_arr[izidx], thetaG_arcmin[j]), out[:,j,:])
#out_std = array(pool.map(std_map, range(1,10001)))
#print out_std.shape
#save(out_dir+'PDFs_Mnu0.%i_z%.1f.npy'%(imnu, z_arr[izidx]), out_std)

pool.close()
print 'DONE DONE DONE'
