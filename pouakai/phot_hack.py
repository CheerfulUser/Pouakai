from core import *
from joblib import Parallel, delayed

cores = 10

files = glob('/home/phys/astronomy/rri38/fli/data/reduced/au_mic/*.gz')

def get_phot(file):
    hdul = fits.open(file)
    wcs = WCS(file)
    cal = cal_photom(data=hdul[0].data,wcs=wcs, header=hdul[0].header,
							 threshold=3,cal_model='ckmodel',
							 brightlim=3,rescale=True)
    name = 'au_mic_phot/' + file.split('/')[-1].split('_')[0] + '_phot.fits'
    hdu = fits.BinTableHDU(data=cal.ap_photom,header=hdul[0].header)
    hdu.writeto(name)

Parallel(cores)(delayed(get_phot)(file) for file in files)

