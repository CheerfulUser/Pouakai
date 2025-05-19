from astropy.io import fits
import pandas as pd
import numpy as np
from glob import glob
from joblib import Parallel, delayed
from tqdm import tqdm

from astropy.coordinates import SkyCoord, get_moon
import astropy.units as u
from astropy.time import Time

import os 
import fnmatch

package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'


moa_darks_dir = '/home/phys/astro8/MJArchive/MOA/DARK/'
moa_flats_dir = '/home/phys/astro8/MJArchive/MOA/FLAT/'
moa_obs_dir = '/home/phys/astro8/MJArchive/MOA/ALERT/'
fli_dir = '/home/phys/astro8/MJArchive/octans/'
fli_dir = '/home/phys/astronomy/zgl12/Pouakai_Test_Data/'

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

def sort_darks(verbose=False,num_core=25, telescope = 'moa'):

	if telescope.lower() == 'moa':
		d = np.append(glob(moa_darks_dir + '*/*.gz'),glob(moa_darks_dir + '*.gz')) # first for year sudbdir, second for file
		dark_files = set(d) 
		try:
			dark_list = pd.read_csv(package_directory + 'cal_lists/moa_dark_list.csv')
		except:
			dark_list = pd.DataFrame(columns=['name','telescope','exptime','jd','date','chip','readodut','filename'])
	else:
		ds = []
		pattern = '*dark*.fit'
		for root, dirs, files in os.walk(fli_dir):
			for file in files:
				if fnmatch.fnmatchcase(file.lower(), pattern.lower()):
					ds.append(os.path.join(root, file))
		dark_files = set(ds)
		try:
			dark_list = pd.read_csv(package_directory + 'cal_lists/bc_dark_list.csv')
		except:
			dark_list = pd.DataFrame(columns=['name','telescope','exptime','jd','date','chip','readodut','filename', 'shape'])

	#dark_list = _kill_old_paths(dark_list,'dark')
	old = set(dark_list['filename'])
	new = dark_files - old
	if verbose: 
		print('Number of new darks: ',len(new))
	files = list(new)
	if len(files) > 0:
		entries = Parallel(num_core)(delayed(dark_info_grab)(file,verbose, telescope) for file in tqdm(files, desc='Processing files'))
		for entry in entries:
			dark_list = dark_list.append(entry, ignore_index=True)
		if telescope.lower() == 'moa':
			dark_list.to_csv(package_directory + 'cal_lists/moa_dark_list.csv',index=False)
		else:
			dark_list.to_csv(package_directory + 'cal_lists/bc_dark_list.csv',index=False)
	if verbose:
		print('Updated darks')

def dark_info_grab(file,verbose = True, telescope = 'moa'):
	entry = {}
	name = file.split('/')[-1].split('.')[0]
	entry['name'] = name

	if telescope.lower() == 'moa':
		try:
			header = fits.open(file)[0].header
			entry['telescope'] = 'MOA-II'
			entry['exptime'] = header['EXPTIME']
			entry['jd'] = (header['JDEND'] + header['JDSTART'])/2
			entry['date'] = header['DATE-OBS'].strip()
			entry['chip'] = header['CHIP']
			entry['readodut'] = -1
		except:
			entry['telescope'] = 'bad'
			entry['exptime'] = 'bad'
			entry['jd'] = 'bad'
			entry['date'] = 'bad'
			entry['chip'] = 'bad'
			entry['readodut'] = 'bad'
	else:
		try:
			header = fits.open(file)[0].header
			entry['telescope'] = 'B&C'
			entry['exptime'] = header['EXPTIME']
			entry['jd'] = header['JD']
			entry['date'] = header['DATE-OBS'].strip()
			entry['chip'] = -1
			entry['readodut'] = header['READOUTM'].replace(" ", "")
			entry['shape'] = fits.open(file)[0].data.shape[0]
		except:
			entry['telescope'] = 'bad'
			entry['exptime'] = 'bad'
			entry['jd'] = 'bad'
			entry['date'] = 'bad'
			entry['chip'] = 'bad'
			entry['readodut'] = 'bad'
			entry['shape'] = 'bad'
	
	entry['filename'] = file
	df = pd.DataFrame([entry])
	return df


def sort_flats(verbose = False, num_core = 25, telescope = 'moa'):
	if telescope.lower() == 'moa':
		f = np.append(glob(moa_flats_dir + '*/*.gz'),glob(moa_flats_dir + '*.gz')) # first for year sudbdir, second for file
		flat_files = set(f)
		try:
			flat_list = pd.read_csv(package_directory + 'cal_lists/moa_flat_list.csv')
		except:
			flat_list = pd.DataFrame(columns=['name','telescope','exptime','jd','date','chip','band','readout','filename','note'])
	else:
		flats = glob(fli_dir + '*/*flat*.fit') # first for day sudbdir, second for file
		Flats = glob(fli_dir + '*/*Flat*.fit') # first for day sudbdir, second for file
		f = np.append(flats,Flats)
		flat_files = set(f) 
		try:
			flat_list = pd.read_csv(package_directory + 'cal_lists/bc_flat_list.csv')
		except:
			flat_list = pd.DataFrame(columns=['name','telescope','exptime','jd','date','chip','band','readout','filename', 'shape' 'note'])
	#flat_list = _kill_old_paths(flat_list,'flat')
	old = set(flat_list['filename'])
	new = flat_files - old

	if verbose: 
		print('Number of new flats: ',len(new))
	files = list(new)
	if len(files) > 0:
		entries = Parallel(n_jobs=num_core)(delayed(flat_info_grab)(file, verbose, telescope) for file in tqdm(files, desc='Processing files'))
		for entry in entries:
			flat_list = flat_list.append(entry, ignore_index=True)
		if telescope.lower() == 'moa':
			flat_list.to_csv(package_directory + 'cal_lists/moa_flat_list.csv',index=False)
		else:
			flat_list.to_csv(package_directory + 'cal_lists/bc_flat_list.csv',index=False)
	if verbose:
		print('Updated flats')


def flat_info_grab(file,verbose=False, telescope = 'moa'):

	name = file.split('/')[-1].split('.')[0]
	entry = {}
	entry['name'] = name

	if telescope.lower() == 'moa':
		hdu = fits.open(file)
		header = hdu[0].header
		data = hdu[0].data
		hdu.close()
		average = np.nanmedian(data)
		if average < 18000:
			note = 'lower'
		elif average > 45000:
			note = 'over'
		else:
			note = 'good'

		entry['band'] = header['COLOUR'].strip()
		entry['chip'] = header['CHIP']
		entry['readout'] = -1
		entry['telescope'] = 'MOA-II'
		entry['exptime'] = header['EXPTIME']
		entry['jd'] = header['JDSTART']
		entry['date'] = header['DATE-OBS'].strip()
		entry['field'] = header['FIELD'].strip()
		entry['filename'] = file
		entry['note'] = note

	else:
		try:
			hdu = fits.open(file)
			header = hdu[0].header
			data = hdu[0].data
			hdu.close()
			average = np.nanmedian(data)
			if average < 18000:
				note = 'lower'
			elif average > 45000:
				note = 'over'
			else:
				note = 'good'

			entry['band'] = header['FILTER'].strip()
			entry['exptime'] = header['EXPTIME']
			entry['telescope'] = 'B&C'
			entry['readout'] = header['READOUTM'].replace(" ", "")
			entry['chip'] = -1
			entry['jd'] = header['JD']
			entry['date'] = header['DATE-OBS'].strip()
			entry['filename'] = file
			entry['note'] = note
			entry['shape'] = data.shape[0]
		except:
			entry['band'] = 'bad'
			entry['chip'] = 'bad'
			entry['telescope'] = 'bad'
			entry['readout'] = 'bad'
			entry['exptime'] = 'bad'
			entry['jd'] = 'bad'
			entry['date'] = 'bad'
			entry['filename'] = file
			entry['note'] = 'bad'
			entry['shape'] = 'bad'

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
		entries = Parallel(num_core)(delayed(dark_info_grab)(file,verbose, telescope) for file in tqdm(files, desc='Processing files'))
		for entry in entries:
			if len(entry) > 0:
				obs_list = obs_list.append(entry, ignore_index=True)
		obs_list.to_csv(package_directory + 'cal_lists/obs_list.csv',index=False)
	if verbose:
		print('Updated obs')

def bias_scrubber(files):
	good = []
	for file in files:
		if 'bias' not in file.lower():
			good += [file]
	return good

def obs_grab_info(file,verbose=False):
	
	entry = {}
	if type(file) == str:
		name = file.split('/')[-1].split('.')[0]
		entry['name'] = name.strip()
		try:
			if telescope.lower() == 'moa':
				header = fits.open(file)[0].header
				entry['telescope'] = 'MOA-II'
				entry['field'] = header['FIELD'].strip()
				entry['chip'] = header['CHIP']
				entry['band'] = header['COLOUR'].strip()
				entry['exptime'] = header['EXPTIME']
				entry['jd'] = (header['JDSTART'] + header['JDEND'])/2
				entry['date'] = header['DATE-OBS'].strip()
				entry['readout'] = -1
			
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
			else:
				header = fits.open(file)[0].header
				entry['telescope'] = 'B&C'
				entry['field'] = header['OBJECT'].replace(' ', '_')
				entry['band'] = header['FILTER'].strip()
				entry['exptime'] = header['EXPTIME']
				entry['jd'] = header['JD']
				entry['date'] = header['DATE-OBS'].strip()
				entry['readout'] = header['READOUTM'].replace(" ", "")
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

		except:
			entry['field'] = 'bad'
			entry['chip'] = 'bad'
			entry['band'] = 'bad'
			entry['readout'] = 'bad'
			entry['telescope'] = 'bad'
			entry['exptime'] = 'bad'
			entry['jd'] = 'bad'
			entry['date'] = 'bad'
			entry['ra'] = 'bad'
			entry['dec'] = 'bad'
			entry['moon_sep'] = 'bad'
			entry['sky'] = 'bad'

		entry['filename'] = file

		df = pd.DataFrame([entry])
		return df

def sort_cals(verbose=True, num_cores = 25, telescope = 'moa'):
	sort_darks(verbose, num_cores, telescope)
	sort_flats(verbose, num_cores, telescope)


if __name__=='__main__':
	sort_darks(verbose=True)
	sort_flats(verbose=True)
	sort_obs(verbose=True)
