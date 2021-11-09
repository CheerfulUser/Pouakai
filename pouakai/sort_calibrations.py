from astropy.io import fits
import pandas as pd
import numpy as np
from glob import glob
from unlzw import unlzw


moa_darks_dir = '/home/phys/astro8/MJArchive/MOA/DARK/'
moa_flats_dir = '/home/phys/astro8/MJArchive/MOA/FLAT/'

def sort_darks():

	dark_list = pd.read_csv('dark_list.csv')
	
	dark_files = set(glob(moa_darks_dir))

	old = set(dark_list['filename'])

	new = dark_files ^ old

	while len(new) > 0:
		entry = {}
		n = new.pop()
		name = n.split('/')[-1].split('.')[0]
		entry['name'] = name
		header = fits.HDUList.fromstring(unlzw(open(n,mode='rb').read()))[0].header
		entry['chip'] = header['CHIP']
		entry['exptime'] = header['EXPTIME']
		entry['jd'] = header['JDSTART']
		entry['date'] = header['DATE-OBS']
		entry['filename'] = n
		
		darktable = darktable.append(entry, ignore_index=True)

	darktable.to_csv('dark_list.csv',index=False)


def sort_darks():

	flat_list = pd.read_csv('flat_list.csv')
	
	flat_files = set(glob(moa_flats_dir))

	old = set(flat_list['filename'])

	new = flat_files ^ old

	while len(new) > 0:
		entry = {}
		n = new.pop()
		name = n.split('/')[-1].split('.')[0]
		entry['name'] = name
		header = fits.HDUList.fromstring(unlzw(open(n,mode='rb').read()))[0].header
		entry['band'] = header['COLOUR']
		entry['chip'] = header['CHIP']
		entry['exptime'] = header['EXPTIME']
		entry['jd'] = header['JDSTART']
		entry['date'] = header['DATE-OBS']
		entry['filename'] = n
		
		darktable = flattable.append(entry, ignore_index=True)

	flattable.to_csv('flat_list.csv',index=False)