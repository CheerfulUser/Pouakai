from astropy.io import fits
import pandas as pd
import numpy as np
from glob import glob


moa_darks_dir = '/home/phys/astro8/rri38/moa/data/DARK/'
moa_flats_dir = '/home/phys/astro8/rri38/moa/data/FLAT/'

def sort_darks(verbose=False):
	dark_files = set(glob(moa_darks_dir + '*.gz'))
	try:
		dark_list = pd.read_csv('cal_lists/dark_list.csv')
		old = set(dark_list['filename'])
		new = dark_files ^ old
	except:
		new = dark_files
		dark_list = pd.DataFrame(columns=['name','chip','exptime','jd','date','filename','note'])

	while len(new) > 0:
		entry = {}
		n = new.pop()
		name = n.split('/')[-1].split('.')[0]
		entry['name'] = name
		try:
			header = fits.open(n)[0].header

			entry['chip'] = header['CHIP']
			entry['exptime'] = header['EXPTIME']
			entry['jd'] = header['JDSTART']
			entry['date'] = header['DATE-OBS']
			entry['filename'] = n
		except:
			print('bad ',n)
		if verbose:
			print('Done ', n)
			print('len n ',len(new))
		
		dark_list = dark_list.append(entry, ignore_index=True)
		dark_list.to_csv('cal_lists/dark_list.csv',index=False)

	


def sort_flats(verbose = False):
	flat_files = set(glob(moa_flats_dir + '*.gz'))
	try:
		flat_list = pd.read_csv('cal_lists/flat_list.csv')	
		old = set(flat_list['filename'])
		new = flat_files ^ old
	except:
		new = flat_files
		flat_list = pd.DataFrame(columns=['name','band','chip','exptime','jd','date','field','filename','note'])

	while len(new) > 0:
		entry = {}
		n = new.pop()
		name = n.split('/')[-1].split('.')[0]
		entry['name'] = name
		header = fits.open(n)[0].header
		data = fits.open(n)[0].data
		average = np.nanmedian(data)
		if average < 18000:
			note = 'lower'
		elif average > 45000:
			note = 'over'
		else:
			note = 'good'

		entry['band'] = header['COLOUR']
		entry['chip'] = header['CHIP']
		entry['exptime'] = header['EXPTIME']
		entry['jd'] = header['JDSTART']
		entry['date'] = header['DATE-OBS']
		entry['field'] = header['FIELD']
		entry['filename'] = n
		entry['note'] = note
		if verbose:
			print('Done ', n)
			print('len n ',len(new))
		
		flat_list = flat_list.append(entry, ignore_index=True)
		flat_list.to_csv('cal_lists/flat_list.csv',index=False)



if __name__=='__main__':
	sort_darks(verbose=True)
	sort_flats(verbose=True)
