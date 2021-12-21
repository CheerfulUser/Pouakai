from astropy.io import fits
import pandas as pd
import numpy as np
from glob import glob


def split_names(files):
	names = [x.split('-')[0] for x in files]
	return names


def make_master_darks(save_location = '/home/users/rri38/data/dark/',verbose=False):
	# make save_location an environment variable
	dark_list = pd.read_csv('cal_lists/dark_list.csv')
	masters = pd.read_csv('cal_lists/master_dark_list.csv')
	names = split_names(dark_list['name'].values)
	all_names = set(names)
	master_names = set(split_names(masters['name'].values))

	new = all_names ^ master_names
	new = list(new)
	new.sort(reverse=True)
	for i in range(len(new)):
		entry = {}
		n = new[i]
		ind = np.array(names) == n
		all_chips = dark_list.iloc[ind]

		for i in range(10):
			i += 1
			chip_ind = all_chips['chip'].values == i
			chip = all_chips.iloc[chip_ind]
			chip_files = chip['filename'].values
			master = []
			for file in chip_files:
				hdu = fits.open(file)[0]
				header = hdu.header
				data = hdu.data
				master += [data]
			master = np.array(master)
			if verbose:
				print('Used ',len(master),' images in median')
			m = np.nanmedian(master,axis=0)
			std = np.nanstd(master,axis=0)
			time = np.nanmean(chip['jd'])
			header['JDSTART'] = time 
			header['MASTER'] = True
			phdu = fits.PrimaryHDU(data = m, header = header)
			ehdu = fits.ImageHDU(data = std, header = header)
			hdul = fits.HDUList([phdu, ehdu])


			letter = file.split('-')[2]
			base_name = file.split('/')[-1].split('.')[0].replace(letter,'m')
			save_name = save_location + base_name + '.fits'
			hdul.writeto(save_name)
		
			entry['name'] = base_name

			entry['chip'] = header['CHIP']
			entry['exptime'] = header['EXPTIME']
			entry['jd'] = time
			entry['date'] = header['DATE-OBS']
			entry['filename'] = save_name
			if len(master) < 3:
				note = 'bad'
			else:
				note = 'good'
			entry['note'] = note
		
			if verbose:
				print('Done ', base_name)
		
			masters = masters.append(entry, ignore_index=True)
			masters.to_csv('cal_lists/master_dark_list.csv',index=False)


def get_master_dark(jd,exptime,chip):
	darks = pd.read_csv('cal_lists/master_dark_list.csv')
	dchips = darks['chip'].values
	dexptime = darks['exptime'].values
	chip_ind = dchips == chip

	exp_ind = dexptime == exptime
	ind = chip_ind & exp_ind
	good = darks.iloc[ind]
	if len(good) > 0:
		djd = good['jd'].values
		diff = jd - djd
		min_ind = np.argmin(abs(diff))
		t_diff = diff[min_ind]
		dark = darks.iloc[min_ind]
		fname = dark['filename'].values
		return fname, t_diff
	else:
		return 'none', -999




def make_master_flats(save_location = '/home/users/rri38/data/flat/',verbose=False):
	# make save_location an environment variable
	flat_list = pd.read_csv('cal_lists/flat_list.csv')
	masters = pd.read_csv('cal_lists/master_flat_list.csv')
	names = split_names(flat_list['name'].values)
	all_names = set(names)
	master_names = set(split_names(masters['name'].values))

	new = all_names ^ master_names
	new = list(new)
	new.sort(reverse=True)
	for i in range(len(new)):
		entry = {}
		n = new[i]
		ind = np.array(names) == n
		all_chips = flat_list.iloc[ind]

		for i in range(10):
			i += 1
			chip_ind = all_chips['chip'].values == i
			chip = all_chips.iloc[chip_ind]
			chip_files = chip['filename'].values
			master = []
			for file in chip_files:
				hdu = fits.open(file)[0]
				header = hdu.header
				data = hdu.data.astype(float)

				saturations = (data > 40000).flatten()
				if sum(saturations) > 100:
					data = data * np.nan

				master += [data]
			master = np.array(master)
			if verbose:
				print('Used ',len(master),' images in median')
			# get dark frame
			fname, tdiff = get_master_dark(chip['jd'].values, chip['exptime'].values, i)
			if verbose:
				print('using dark frame ',fname)
			try:
				dark = fits.open(fname)[0].data
				master = master - dark
			except:
				m = '!!! Warning: No dark found !!!'
				print(m)
			m = np.nanmedian(master,axis=0)
			std = np.nanstd(master,axis=0)
			time = np.nanmean(chip['jd'])
			header['JDSTART'] = time 
			header['MASTER'] = True
			phdu = fits.PrimaryHDU(data = m, header = header)
			ehdu = fits.ImageHDU(data = std, header = header)
			hdul = fits.HDUList([phdu, ehdu])


			letter = file.split('-')[2]
			base_name = file.split('/')[-1].split('.')[0].replace(letter,'m')
			save_name = save_location + base_name + '.fits'
			hdul.writeto(save_name)
		
			entry['name'] = base_name

			entry['chip'] = header['CHIP']
			entry['exptime'] = header['EXPTIME']
			entry['jd'] = time
			entry['date'] = header['DATE-OBS']
			entry['filename'] = save_name
			entry['dark_file'] = fname 
			entry['time_diff'] = tdiff
			if (np.nanmedian(m) < 15000):
				note = 'bad'
			else:
				note = 'good'
			entry['note'] = note
			
			field = header['FIELD']
			if 'flat_round' in field:
				flat_type = 'dome'
			else:
				flat_type = 'sky'
			
			entry['flat_type'] = flat_type

			if verbose:
				print('Done ', base_name)
		
			masters = masters.append(entry, ignore_index=True)
			masters.to_csv('cal_lists/master_flat_list.csv',index=False)


