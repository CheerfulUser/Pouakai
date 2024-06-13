from astropy.io import fits
import pandas as pd
import numpy as np
from glob import glob
from joblib import Parallel, delayed

from astropy.coordinates import SkyCoord, get_moon
import astropy.units as u
from astropy.time import Time

fli_dir = '/home/phys/astro8/MJArchive/octans/'

def sort_darks(verbose=False,num_core=25):
	darks = glob(fli_dir + '*/*dark*.fit') 
	Darks = np.append(glob(fli_dir + '*/*Dark*.fit'), glob(fli_dir + '*/*DARK*.fit'))
	d = np.append(darks,Darks)
	dark_files = set(d) 
	dark_list = pd.read_csv('cal_lists/dark_list.csv')
	old = set(dark_list['filename'])
	new = dark_files - old
	if verbose: 
		print('Number of new darks: ',len(new))
	files = list(new)
	if len(files) > 0:
		entries = Parallel(num_core)(delayed(dark_info_grab)(file,verbose) for file in files)
		for entry in entries:
			dark_list = pd.concat([dark_list, entry], ignore_index=True)
		dark_list.to_csv('cal_lists/dark_list.csv',index=False)
	if verbose:
		print('Updated darks')


def dark_info_grab(file,verbose = True):
	entry = {}
	name = file.split('/')[-1].split('.')[0]
	entry['name'] = name
	try:
		header = fits.open(file)[0].header
		t = header['TELESCOP']
		if len(t) > 0:
			entry['telescope'] = header['TELESCOP'].strip()
		else:
			entry['telescope'] = 'B&C'
		entry['exptime'] = header['EXPTIME']
		entry['jd'] = header['JD']
		entry['date'] = header['DATE-OBS'].strip()
	except:
		print('!!! bad ',file)
		entry['telescope'] = 'bad'
		entry['exptime'] = 'bad'
		entry['jd'] = 'bad'
		entry['date'] = 'bad'
	
	entry['filename'] = file
	if verbose:
		print('Done ', file)
	df = pd.DataFrame([entry])
	return df

	


def sort_flats(verbose = False, num_core = 25):
	flats = glob(fli_dir + '*/*flat*.fit') # first for day sudbdir, second for file
	Flats = glob(fli_dir + '*/*Flat*.fit') # first for day sudbdir, second for file
	f = np.append(flats,Flats)
	flat_files = set(f) 
	flat_list = pd.read_csv('cal_lists/flat_list.csv')	
	old = set(flat_list['filename'])
	new = flat_files - old

	if verbose: 
		print('Number of new flats: ',len(new))
	files = list(new)
	if len(files) > 0:
		entries = Parallel(num_core)(delayed(flat_info_grab)(file,verbose) for file in files)
		for entry in entries:
			flat_list = pd.concat([flat_list, entry], ignore_index=True)
		flat_list.to_csv('cal_lists/flat_list.csv',index=False)
	if verbose:
		print('Updated flats')


def flat_info_grab(file,verbose=False):

	name = file.split('/')[-1].split('.')[0]
	entry = {}
	entry['name'] = name
	try:
		hdu = fits.open(file)[0]
		header = hdu.header
		data = hdu.data
		average = np.nanmedian(data)
		if average < 18000:
			note = 'lower'
		elif average > 45000:
			note = 'over'
		else:
			note = 'good'

		
	
		entry['band'] = header['FILTER'].strip()
		t = header['TELESCOP']
		if len(t) > 0:
			entry['telescope'] = header['TELESCOP'].strip()
		else:
			entry['telescope'] = 'B&C'
		entry['exptime'] = header['EXPTIME']
		entry['jd'] = header['JD']
		entry['date'] = header['DATE-OBS'].strip()
		entry['filename'] = file
		entry['note'] = note
	except:
		entry['band'] = 'bad'
		entry['telescope'] = 'bad'
		entry['exptime'] = -999
		entry['jd'] = -999
		entry['date'] = 'bad'
		entry['filename'] = 'bad'
		entry['note'] = 'bad'
	
	if verbose:
		print('Done ', name)
	df = pd.DataFrame([entry])
	return df

def bias_scrubber(files):
	good = []
	for file in files:
		if 'bias' not in file.lower():
			good += [file]
	return good


def sort_obs(verbose=False,num_core = 25):
	all_ims = glob(fli_dir + '*/*.fit')
	a = set(all_ims)
	flats = glob(fli_dir + '*/*flat*.fit')
	Flats = glob(fli_dir + '*/*Flat*.fit')
	f = np.append(flats,Flats)

	darks = glob(fli_dir + '*/*dark*.fit')
	Darks = np.append(glob(fli_dir + '*/*Dark*.fit'), glob(fli_dir + '*/*DARK*.fit'))
	d = np.append(darks,Darks)
	cal = np.append(d,f)
	cals = set(cal) 

	obs_files = a - cals

	obs_list = pd.read_csv('cal_lists/obs_list.csv')
	old = set(obs_list['filename'].values)
	new = obs_files ^ old
	files = list(new)
	files = bias_scrubber(files)
	if verbose: 
		print('Number of new obs: ',len(files))

	if len(files) > 0:
		entries = Parallel(num_core)(delayed(obs_grab_info)(file,verbose) for file in files)
		for entry in entries:
			if len(entry) > 0:
				obs_list = obs_list.append(entry, ignore_index=True)
		obs_list.to_csv('cal_lists/obs_list.csv',index=False)
	if verbose:
		print('Updated obs')


def obs_grab_info(file,verbose=False):
	
	entry = {}
	if type(file) == str:
		name = file.split('/')[-1].split('.')[0]
		entry['name'] = name.strip()
		try:
			header = fits.open(file)[0].header
			entry['object'] = header['OBJECT'].strip()
			t = header['TELESCOP']
			if len(t) > 0:
				entry['telescope'] = header['TELESCOP'].strip()
			else:
				entry['telescope'] = 'B&C'
			try:
				entry['band'] = header['FILTER'].strip()
			except:
				print('!!! Filter not defined !!!')
				entry['band'] = 'bad'
			entry['exptime'] = header['EXPTIME']
			entry['jd'] = header['JD']
			entry['date'] = header['DATE-OBS'].strip()
			try:
				ra = header['RA'].strip()
				dec = header['DEC'].strip()
				c = SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
				t = Time(header['JD'],format='jd')
				moon = get_moon(t)
				sep = moon.separation(c)
				entry['ra'] = ra
				entry['dec'] = dec
				entry['moon_sep'] = sep.deg
			except:
				print('No coordinates defined')

			entry['sky'] = np.nanmedian(fits.open(file)[0].data)
			entry['image_type'] = header['IMAGETYP']

			entry['filename'] = file
		except:
			entry['object'] = 'bad'
			entry['telescope'] = 'bad'
			entry['band'] = 'bad'
			entry['exptime'] = 'bad'
			entry['jd'] = 'bad'
			entry['date'] = 'bad'
			entry['ra'] = 'bad'
			entry['dec'] = 'bad'
			entry['moon_sep'] = 'bad'
			entry['sky'] = 'bad'
			entry['image_type'] = 'bad'
			entry['filename'] = file
		
		if verbose:
			print('Done ', name)

		df = pd.DataFrame([entry])
		return df

def sort_cals(verbose=True):
	sort_darks(verbose)
	sort_flats(verbose)

if __name__=='__main__':
	sort_darks(verbose=True)
	sort_flats(verbose=True)
	sort_obs(verbose=True)
