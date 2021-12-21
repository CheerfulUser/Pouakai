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
	new = list(nmew)
	for i in range(len(new)):
		entry = {}
		n = new[i]
		ind = np.array(names) == n
		all_chips = dark_list.iloc[ind]

		for i in range(10):
			i += 1
			chip_ind = all_chips['chip'].value == i
			chip = all_chips.iloc[chip_ind]
			chip_files = chip['filename'].values
			master = []
			for file in chip_files:
				hdu = fits.open(file)[0]
				header = hdu.header
				data = hdu.data
				master += [data]
			master = np.array(master)
			m = np.nanmedian(master,axis=0)
			std = np.nansrd(master,axis=0)
			time = np.nanmean(chip['jd'])
			header['JDSTART'] = time 
			header['MASTER'] = True
			phdu = fits.PrimaryHDU(data = m, header = header)
			ehdu = fits.ImageHDU(data = std, header = header)
			hdul = fits.HDUList([phdu, ehdu])


			letter = file.split('-')[2]
			base_name = file.replace(letter,'m')
			save_name = save_location + base_name + 'fits.gz'
			hdul.writeto(save_name)
		
			entry['name'] = base_name

			entry['chip'] = header['CHIP']
			entry['exptime'] = header['EXPTIME']
			entry['jd'] = time
			entry['date'] = header['DATE-OBS']
			entry['filename'] = save_name
		
			if verbose:
				print('Done ', n)
		
			masters = masters.append(entry, ignore_index=True)
			masters.to_csv('cal_lists/master_dark_list.csv',index=False)




