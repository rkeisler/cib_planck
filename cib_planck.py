import numpy as np

datadir = './data/'

def get_cl_cibcib(l_interp):
    '''
    Returns the temperature spectra due to the CIB (clustered+shot)
    for different Planck HFI bands.  Emailed from G. Lagache to RK.
    More info:

    The measurements: there is one file per frequency.
    Note that the first two low-ell points have to be considered 
    as upper limits due to Galactic dust residuals.

    The best fit extended HOD model. 

    They are Cells in Jy^2/sr. To convert them in microK_CMB^2, 
    you just need to divide the Cell by the square of the 
    conversion factor given in Table 6 of arXiv:1303.5070.

    The shot noise is not included in the model file.
    '''
    from numpy import interp
    from pandas.io.parsers import read_csv
    from itertools import combinations_with_replacement
    d=read_csv(datadir+'all_best_fit.txt',sep=' ')
    # The raw data are in Jy^2/sr.
    # We convert them to uK^2 below.
    bands=['143','217','353','545','857']
    l_csv = d['ell']
    cl_shot = get_cl_shot()
    output={}
    for band1,band2 in combinations_with_replacement(bands, 2):
        this_key = 'Cl'+band1+'x'+band2
        if not(this_key in d.columns): 
            this_key = 'Cl'+band2+'x'+band1
        this_cl=d[this_key]
        this_cl_interp = interp(l_interp, l_csv, this_cl)
        this_cl_interp /= jy_per_sr_PER_uk_cmb(band1)
        this_cl_interp /= jy_per_sr_PER_uk_cmb(band2)
        this_cl_interp += cl_shot[band1,band2]
        output[band1,band2]=this_cl_interp
        output[band2,band1]=this_cl_interp
    return output 


def get_cl_shot():
    '''
    Returns the shot power C^{shot}, from CIB+radio, 
    for Planck HFI bands.
    '''
    # from Tables 6 & 7 of 1309.0382
    # We'll enter them in Jy**2/sr, 
    # then convert to (uK_CMB)**2 at the end.
    cl={}
    cl['857','857'] = 5628 + 4.28
    cl['857','545'] = 2655 + 2.28
    cl['857','353'] = 913 + 2.10
    cl['857','217'] = 216 + 1.53
    cl['857','143'] = 56 + 2.38
    
    cl['545','545'] = 1454 + 2.86
    cl['545','353'] = 543 + 2.59
    cl['545','217'] = 135 + 1.92
    cl['545','143'] = 35 + 2.86

    cl['353','353'] = 225 + 3.28
    cl['353','217'] = 59 + 2.4
    cl['353','143'] = 15 + 3.57

    cl['217','217'] = 16 + 3.12
    cl['217','143'] = 4.3 + 3.68

    cl['143','143'] = 1.2 + 6.05

    # symmetrize the dictionary
    cl2={}
    for k,v in cl.iteritems():
        cl2[k[0],k[1]]=v
        cl2[k[1],k[0]]=v
    cl=cl2
        
    # convert to (uK_CMB**2)
    for k,v in cl.iteritems():
        cl[k] = v/jy_per_sr_PER_uk_cmb(k[0])/jy_per_sr_PER_uk_cmb(k[1])
    return cl


def get_cl_phicib(l_interp):
    '''
    Returns the cross-spectrum between the CMB lensing potential phi
    and the CIB temperature anisotropies, C^{\phi T}.  The files 
    were emailed from G. Lagache to RK and are "lensingxCIB prediction 
    for the extended halo model".
    '''
    bands=['143','217','353','545','857']
    output={}
    for band in bands:
        tmp=np.loadtxt(datadir+'lensing_HOD/lensing_'+band+'.txt')
        l_tmp = tmp[:,0]
        cl_tmp = tmp[:,1]
        cl_tmp /= jy_per_sr_PER_uk_cmb(band)
        output[band] = np.interp(l_interp, l_tmp, cl_tmp) / ((l_interp+1.)**3.)
    return output


def jy_per_sr_PER_uk_cmb(band):
    '''
    A conversion utility to go from (Jy per sr) to uK-CMB
    for the Planck HFI bands.  From Table 6 1303.5070
    '''
    # Given in MJy/sr per K_CMB, which is the same as 
    # Jy/sr per uK_CMB.
    xdust = {'143':358.04, '217':415.465, 
             '353':246.543, '545':49.59, '857':2.09}
    return xdust[band]

