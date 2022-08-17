from astropy.io import fits
import pandas as pd
import numpy as np
from glob import glob
from astropy.coordinates import SkyCoord, get_moon
import astropy.units as u
from astropy.time import Time


moa_obs_dir = '/home/phys/astro8/MJArchive/MOA/ALERT/'

def sort_obs(verbose=False):
	obs_files = set(glob(moa_obs_dir + '*.gz'))
	#try:
	obs_list = pd.read_csv('cal_lists/obs_list.csv')
	old = set(obs_list['filename'])
	new = obs_files ^ old
	print(len(new))

	while len(new) > 0:
		entry = {}
		n = new.pop()
		print(n)
		name = n.split('/')[-1].split('.')[0]
		entry['name'] = name
		try:
			header = fits.open(n)[0].header

			entry['chip'] = header['CHIP']
			entry['exptime'] = header['EXPTIME']
			entry['jd'] = header['JDSTART']
			entry['date'] = header['DATE-OBS']
			
			ra = header['RA'].strip()
			dec = header['DEC'].strip()
			c = SkyCoord(ra,dec,units=(u.hourangle,u.deg))
			t = Time(header['JDSTART'],format='jd')
			moon = get_moon(t)
			sep = moon.separation(c)
			entry['ra'] = ra
			entry['dec'] = dec
			entry['moon_sep'] = sep
			entry['sky'] = np.nanmedian(fits.open(n)[0].data)

			entry['filename'] = n

		except:
			print('bad ',n)
		if verbose:
			print('Done ', n)
			print('len n ',len(new))
		
		obs_list = obs_list.append(entry, ignore_index=True)
		obs_list.to_csv('cal_lists/obs_list.csv',index=False)

if __name__=='__main__':
	sort_obs(verbose=True)