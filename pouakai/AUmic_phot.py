from aperture_photom import cal_photom
import numpy as np
from glob import glob
from astropy.io import fits
from astropy.wcs import WCS

files = glob('/home/phys/astronomy/rri38/fli/data/reduced/cal/*.gz')


def run_phot(file):
    base_name = file.split('/')[-1].split('_')[0]
    hdul = fits.open(file)
    wcs = WCS(file)
    header = hdul[0].header
    cal = cal_photom(data=hdul[0].data,wcs=wcs, header=header,
                                threshold=3,cal_model='ckmodel',ax=None,
                                brightlim=10,rescale=False)

        

    name = '/home/phys/astronomy/rri38/fli/data/reduced/phot_table' + base_name + '_phot.fits'
    ind = cal.ap_photom.mag.values < cal.maglim3
    tab = cal.ap_photom.iloc[ind]
    ra,dec = wcs.all_pix2world(tab['xcenter'].values,tab['ycenter'].values,0)
    tab['ra'] = ra
    tab['dec'] = dec
    print(list(tab.keys()))
    rec = np.rec.array([tab['gaiaID'].values,tab['xcenter'].values,tab['ycenter'].values,tab['ra'].values,tab['dec'].values,tab['counts'].values,tab['e_counts'].values,tab['mag'].values,tab['e_mag'].values,tab['snr'].values],
                        formats='int,float32,float32,float32,float32,float32,float32,float32,float32,float32',
                        names='gaiaID,xcenter,ycenter,ra,dec,counts,counts_e,mag,mag_e,snr' )

    hdu = fits.BinTableHDU(data=rec,header=header)

    hdu.writeto(name)

