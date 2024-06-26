from astropy.io import fits
import pandas as pd
import numpy as np
from glob import glob
from joblib import Parallel, delayed

from astropy.coordinates import SkyCoord, get_moon
import astropy.units as u
from astropy.time import Time

import os 
package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'


moa_darks_dir = '/home/phys/astro8/MJArchive/MOA/DARK/'
moa_flats_dir = '/home/phys/astro8/MJArchive/MOA/FLAT/'
moa_obs_dir = '/home/phys/astro8/MJArchive/MOA/ALERT/'

def _update_paths(files,csv,cal_type):
	csv = csv.drop_duplicates(subset=['name'])
	for file in files:
		name = file.split('/')[-1].split('.')[0]
		ind = np.where(name == csv.name.values)
		if csv.iloc[ind]['filename'].values != file:
			print('!!!! Changing path ' + name)
			csv.iloc[ind]['filename'] = file
	csv.to_csv(package_directory + 'cal_lists/{}_list.csv'.format(cal_type),index=False)
	return csv
def _kill_old_paths(csv,cal_type):
	ind = []
	for i in range(len(csv)):
		try:
			hdu = fits.open(csv.iloc[i]['filename'])
			ind += [i]
		except:
			print('bad path')
	csv_new = csv.iloc[ind]
	csv_new.to_csv(package_directory + 'cal_lists/{}_list.csv'.format(cal_type),index=False)
	return csv_new

def sort_darks(verbose=False,num_core=25):
	f = np.append(glob(moa_darks_dir + '*/*.gz'),glob(moa_darks_dir + '*.gz')) # first for year sudbdir, second for file

	dark_files = set(f) 
	dark_list = pd.read_csv(package_directory + 'cal_lists/dark_list.csv')
	#dark_list = _kill_old_paths(dark_list,'dark')
	old = set(dark_list['filename'])
	new = dark_files - old
	if verbose: 
		print('Number of new darks: ',len(new))
	files = list(new)
	if len(files) > 0:
		entries = Parallel(num_core)(delayed(dark_info_grab)(file,verbose) for file in files)
		for entry in entries:
			dark_list = dark_list.append(entry, ignore_index=True)
		dark_list.to_csv(package_directory + 'cal_lists/dark_list.csv',index=False)
	if verbose:
		print('Updated darks')


def dark_info_grab(file,verbose = True):
	entry = {}
	name = file.split('/')[-1].split('.')[0]
	entry['name'] = name
	try:
		header = fits.open(file)[0].header

		entry['chip'] = header['CHIP']
		entry['exptime'] = header['EXPTIME']
		entry['jd'] = header['JDSTART']
		entry['date'] = header['DATE-OBS'].strip()
	except:
		print('!!! bad ',file)
		entry['chip'] = 'bad'
		entry['exptime'] = 'bad'
		entry['jd'] = 'bad'
		entry['date'] = 'bad'
	
	entry['filename'] = file
	if verbose:
		print('Done ', file)
	df = pd.DataFrame([entry])
	return df

	


def sort_flats(verbose = False, num_core = 25):

	f = np.append(glob(moa_flats_dir + '*/*.gz'),glob(moa_flats_dir + '*.gz')) # first for year sudbdir, second for file
	flat_files = set(f)

	flat_list = pd.read_csv(package_directory + 'cal_lists/flat_list.csv')	
	#flat_list = _kill_old_paths(flat_list,'flat')
	old = set(flat_list['filename'])
	new = flat_files - old

	if verbose: 
		print('Number of new flats: ',len(new))
	files = list(new)
	if len(files) > 0:
		entries = Parallel(num_core)(delayed(flat_info_grab)(file,verbose) for file in files)
		for entry in entries:
			flat_list = flat_list.append(entry, ignore_index=True)
		flat_list.to_csv(package_directory + 'cal_lists/flat_list.csv',index=False)
	if verbose:
		print('Updated flats')


def flat_info_grab(file,verbose=False):

	name = file.split('/')[-1].split('.')[0]
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

	entry = {}
	entry['name'] = name
	entry['band'] = header['COLOUR'].strip()
	entry['chip'] = header['CHIP']
	entry['exptime'] = header['EXPTIME']
	entry['jd'] = header['JDSTART']
	entry['date'] = header['DATE-OBS'].strip()
	entry['field'] = header['FIELD'].strip()
	entry['filename'] = file
	entry['note'] = note
	if verbose:
		print('Done ', name)
	df = pd.DataFrame([entry])
	return df

def sort_obs(verbose=False,num_core = 25):
	f = np.append(glob(moa_obs_dir + '*/*.gz'),glob(moa_obs_dir + '*.gz')) # first for year sudbdir, second for file
	obs_files = set(f)
	#try:
	obs_list = pd.read_csv(package_directory + 'cal_lists/obs_list.csv')
	old = set(obs_list['filename'].values)
	new = obs_files - old
	if verbose: 
		print('Number of new obs: ',len(new))
	files = list(new)
	if len(files) > 0:
		entries = Parallel(num_core)(delayed(dark_info_grab)(file,verbose) for file in files)
		for entry in entries:
			if len(entry) > 0:
				obs_list = obs_list.append(entry, ignore_index=True)
		obs_list.to_csv(package_directory + 'cal_lists/obs_list.csv',index=False)
	if verbose:
		print('Updated obs')


def obs_grab_info(file,verbose=False):
	
	entry = {}

	if type(file) == str:
		name = file.split('/')[-1].split('.')[0]
		entry['name'] = name.strip()
		#try:
		header = fits.open(file)[0].header
		print(header['FIELD'])
		entry['field'] = header['FIELD'].strip()
		entry['chip'] = header['CHIP']
		entry['band'] = header['COLOUR']
		entry['exptime'] = header['EXPTIME']
		entry['jd'] = header['JDSTART']
		entry['date'] = header['DATE-OBS'].strip()
		
		ra = header['RA'].strip()
		dec = header['DEC'].strip()
		c = SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
		t = Time(header['JDSTART'],format='jd')
		moon = get_moon(t)
		sep = moon.separation(c)
		entry['ra'] = ra
		entry['dec'] = dec
		entry['moon_sep'] = sep.deg
		entry['sky'] = np.nanmedian(fits.open(file)[0].data)

		entry['filename'] = file
		'''except:
			entry['field'] = 'bad'
			entry['chip'] = 'bad'
			entry['exptime'] = 'bad'
			entry['jd'] = 'bad'
			entry['date'] = 'bad'
			entry['ra'] = 'bad'
			entry['dec'] = 'bad'
			entry['moon_sep'] = 'bad'
			entry['sky'] = 'bad'

			entry['filename'] = file'''
		
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
