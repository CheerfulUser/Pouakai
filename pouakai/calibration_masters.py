from astropy.io import fits
from astropy.time import Time
import pandas as pd
import numpy as np
from glob import glob
import os
from copy import deepcopy
from joblib import Parallel, delayed
from tqdm import tqdm
import gc
from hidden_prints import HiddenPrints

import matplotlib.pyplot as plt

import time

import warnings
warnings.filterwarnings("ignore")

package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'

def split_names(files):
	names = [str(x).split('-')[0] for x in files]
	return names

def assign_master_name(darks):
	master_name = []
	darks = darks.reset_index(drop=True)
	for i in range(len(darks)):
		dark = darks.iloc[i]
		master_name += [f"dark_{dark['jd']}_{dark['exptime']}s"]
	
	return master_name

def make_master_darks(save_location = '/home/phys/astronomy/rri38/moa/data/master/dark/', time_frame = 1, num_cores=25,verbose=False, telescope = 'moa', dark_tolerance = 1):
	if telescope.lower() == 'moa':
		save_location = '/home/phys/astronomy/rri38/moa/data/master/dark/'
		dark_list = pd.read_csv(package_directory + 'cal_lists/moa_dark_list.csv')
		ind = dark_list['chip'] != 'bad'
		try:
			masters = pd.read_csv(package_directory + 'cal_lists/moa_master_dark_list.csv')
		except:
			masters = pd.DataFrame(columns=['name','telescope','exptime','jd','date','chip','readout','filename', 'nimages', 'note'])
		dark_list = dark_list[ind]
	else:
		save_location = '/home/phys/astronomy/rri38/fli/data/master/dark/'
		dark_list = pd.read_csv(package_directory + 'cal_lists/bc_dark_list.csv')
		try:
			masters = pd.read_csv(package_directory + 'cal_lists/bc_master_dark_list.csv')
		except:
			masters = pd.DataFrame(columns=['name','telescope','exptime','jd','date','chip','readout','filename', 'nimages', 'shape', 'note'])

	if telescope.lower() == 'bc':
		names = assign_master_name(dark_list)
		all_names = set(names)
		master_names = set(masters['name'].values)
	else:
		names = split_names(dark_list['name'].values)
		all_names = set(names)
		master_names = set(split_names(masters['name'].values))

	new = all_names - master_names
	new = list(new)
	new.sort(reverse=True)
	indexer = np.arange(len(new))
	print('Number of new dark entries:', len(new))
	print(save_location)
	if len(new) > 0:
		if num_cores > 1:
			entries = Parallel(n_jobs=num_cores)(delayed(dark_processing)(index,names, new,dark_list,time_frame,save_location,verbose,telescope, dark_tolerance)  for index in tqdm(indexer, desc='Processing files'))
			for entry in entries:
				if entry is not None:
					masters = pd.concat([masters,entry],ignore_index=True)
		else:
			entries = []
			for index in tqdm(indexer, desc='Processing files'):
				entry = dark_processing(index,names, new,dark_list,time_frame,save_location,verbose,telescope, dark_tolerance)
				masters = pd.concat([masters,entry],ignore_index=True)
	
	masters = masters.reset_index(drop=True)

	if telescope.lower() == 'moa':
		masters.to_csv(package_directory + 'cal_lists/moa_master_dark_list.csv',index=False)
	else:
		masters.to_csv(package_directory + 'cal_lists/bc_master_dark_list.csv',index=False)

def _bc_darks(files, save_location):
	master = []
	entry = {}
	try:
		# print('Initial')
		jds = []
		exps = []
		for file in files:
			# print('Files')
			hdu = fits.open(file)[0]
			header = hdu.header
			data = hdu.data
			master.append(data)
			jds.append(float(header['JD']))
			exps.append(float(header['EXPTIME']))
		# print('Appended')
		master = np.array(master)
		m = np.nanmedian(master,axis=0)
		std = np.nanstd(master,axis=0)
		time = np.nanmean(jds)
		header['JDSTART'] = time
		# print('Header')

		utc = Time(time, format='jd').utc.strftime('%Y%m%d')
		header['MASTER'] = True
		phdu = fits.PrimaryHDU(data = m, header = header)
		ehdu = fits.ImageHDU(data = std, header = header)
		hdul = fits.HDUList([phdu, ehdu])
		# print('HDU')

		name = 'master_dark_' + str(utc) + '_' + str(int(np.nanmedian(exps)))

		save_name = save_location + name + '.fits'
		hdul.writeto(save_name,overwrite=True)
		# print('Save')
		compress = 'gzip -f ' + save_name
		with HiddenPrints():
			os.system(compress)
		# print('Compress')
		entry['name'] = name
		entry['telescope'] = 'B&C'
		entry['readout'] = header['READOUTM'].replace(" ", "")
		entry['exptime'] = header['EXPTIME']
		entry['chip'] = -1
		entry['jd'] = time
		entry['date'] = header['DATE-OBS']
		entry['nimages'] = len(master)
		entry['filename'] = save_name + '.gz'
		entry['shape'] = m.shape[0]
		# print('Entry')
		if len(master) < 3:
			# print('BC Darks bad')
			note = 'bad'
		else:
			note = 'good'
		entry['note'] = note
		entry = pd.DataFrame([entry])
			
		return entry
	except:
		# raise Exception('Something went wrong...')
		print('something went wrong...')


def dark_processing(index, names, new, dark_list,time_frame,save_location,verbose,telescope, dark_tolerance):
    
	entries = {}
	name = new[index]
	ind = np.array(names) == name
	all_chips = dark_list.iloc[ind]
 
 
	valid_jd = all_chips.jd.values
	valid_jd = valid_jd[valid_jd != 'bad']  # Remove 'bad' entries

	try:
		t = np.int64(valid_jd.astype(float))[0]
		exptime = int(all_chips.exptime.values.astype(float)[0])
		times = pd.to_numeric(dark_list['jd'].values, errors='coerce')
	except:
		return None

	tind = (t - times >= -time_frame) & (t - times <= time_frame)
	
 
	converted_exptime = pd.to_numeric(dark_list['exptime'].values, errors='coerce')
	expind = converted_exptime == exptime

	if expind.all() == False:
		expind = abs(converted_exptime - exptime) <= dark_tolerance

	if telescope.lower() == 'bc':
		shapes = all_chips['shape'].values.astype(int)[0]
		converted_shapes = pd.to_numeric(dark_list['shape'].values, errors='coerce')
		shapeind = converted_shapes == shapes
		ind = tind & expind & shapeind
	else:
		ind = tind & expind

	files = dark_list['filename'].values[ind]
	jds = times[ind]

	if telescope.lower() == 'moa':
		return _moa_darks(all_chips, save_location, verbose)
	else:
		return _bc_darks(files, save_location)

def _moa_darks(all_chips, save_location, verbose):
	# name,telescope,exptime,jd,date,chip,readodut,filename,nimages,note
	entries = pd.DataFrame(columns = ['name','telescope','exptime','jd','date','chip','readodut','filename','nimages','note'])
	for j in range(10):
		try:
			j += 1
			entry = {}
			chip_ind = all_chips['chip'].values.astype(int) == j
			chip = all_chips.iloc[chip_ind]
			chip_files = chip['filename'].values
			master = []
			for file in chip_files:
				with fits.open(file) as hdu:
					header = hdu[0].header
					data = hdu[0].data
				master.append(data)
			master = np.array(master)
			m = np.nanmedian(master,axis=0)
			std = np.nanstd(master,axis=0)
			time = np.nanmean(chip['jd'].astype(float))

			header['JD'] = time 
			header['MASTER'] = True
			phdu = fits.PrimaryHDU(data = m, header = header)
			ehdu = fits.ImageHDU(data = std, header = header)
			hdul = fits.HDUList([phdu, ehdu])

			letter = file.split('-')[2]
			base_name = file.split('/')[-1].split('.')[0].replace(letter,'m')
			save_name = save_location + base_name + '.fits'
			hdul.writeto(save_name,overwrite=True)
			with HiddenPrints():
				compress = 'gzip -f ' + save_name
				os.system(compress)

			if len(master) < 3:
				print('MOA Darks bad')
				note = 'bad'
			else:
				note = 'good'
			
			entry = [base_name, 'MOA-II', header['EXPTIME'], time, header['DATE-OBS'], header['CHIP'], -1, save_name + '.gz', len(master), note]

			if len(entries) == 0:
				entries = pd.DataFrame([entry], columns = ['name','telescope','exptime','jd','date','chip','readodut','filename','nimages','note'])
			else:
				entries = pd.concat([entries, pd.DataFrame([entry], columns = ['name','telescope','exptime','jd','date','chip','readodut','filename','nimages','note'])])
			gc.collect()
			del data
			del master
		except:
			pass
		
	return entries

def get_master_dark(jd,exptime, chip, readout, strict=True,exp_tol=1, telescope = 'moa', shape = 2048):
	"""
	ytdhgvj
	"""
	
 
	if telescope.lower() == 'moa':
		try:
			darks = pd.read_csv(package_directory + 'cal_lists/moa_master_dark_list.csv')
		except:
			darks = pd.DataFrame(columns=['name','telescope','exptime','jd','date','chip','readodut','filename', 'nimages', 'note'])
	else:
		try:
			darks = pd.read_csv(package_directory + 'cal_lists/bc_master_dark_list.csv')
		except:
			darks = pd.DataFrame(columns=['name','telescope','exptime','jd','date','chip','readodut','filename', 'nimages', 'shape', 'note'])
	if strict:
		ind = darks['note'].values.astype(str) == 'good'
		darks = darks.iloc[ind]

	dchip = darks['chip'].values
	dreadout = darks['readout'].values
	chip_ind = dchip.astype(int) == int(chip)
	readout_ind = dreadout == readout
	
	dexptime = darks['exptime'].values
	exp_ind = abs(np.int64(dexptime.astype(float)) - int(float(exptime))) < exp_tol
	# print('ZZZZ')
 
	if telescope.lower() == 'bc':
		dshape = darks['shape'].values
		shape_ind = dshape.astype(int) == int(shape)
		good = darks.iloc[exp_ind & chip_ind & readout_ind & shape_ind]
		# print('EXP. IND.', np.nansum(exp_ind))
	else:
		good = darks.iloc[exp_ind & chip_ind & readout_ind]
  
	if len(good) > 0:
		djd = good['jd'].values
		djd = djd.astype(float)
		diff = jd - djd
		diff = diff.astype(float)
		min_ind = np.argmin(abs(diff))
		t_diff = diff[min_ind]
		dark = good.iloc[min_ind]
		fname = dark['filename']#.values[0]

		if abs(t_diff) < exp_tol:
			return fname, t_diff
		else:
			return 'none', -999	
	else:
		return 'none', -999

def cut_bad_reductions(table):
	"""
	remove all reductions with no good dark frames to see if anything better can be done
	"""
	ind = table['time_diff'].values == -999
	bad_names = set(split_names(table['name'].iloc[ind]))
	names = split_names(table['name'])
	bad_names = list(bad_names)
	tab = deepcopy(table)
	for i in range(len(bad_names)):
		print('Dropping ' + bad_names[i])
		inds = np.where(names != bad_names[i])[0]
		tab = tab.iloc[inds]
		names = split_names(tab['name'])
	return tab

def _bc_master_flats(save_location, time_frame, num_cores, verbose, exp_tol = 1):
	flat_list = pd.read_csv(package_directory + 'cal_lists/bc_flat_list.csv')
	try:
		masters = pd.read_csv(package_directory + 'cal_lists/bc_master_flat_list.csv')
	except:
		masters = pd.DataFrame(columns=['name','telescope','exptime','jd','date','band','chip','readout','filename','nimages','flattype', 'shape', 'note'])
	
	ind = (flat_list['note'] == 'good')
	flat_list = flat_list[ind]
	
	flat_list['band'] = flat_list['band'].str.strip()
	flat_list['readout'] = flat_list['readout'].str.strip()
	readouts = flat_list['readout'].values
	times = np.int64(np.round(flat_list['jd'].values.astype(float)))
	names = []
	for i in range(len(times)):
		names += [('flat_' + times.astype(str)[i] + '_' + str(time_frame) + 'd_' + flat_list['band'].values[i]) + '_' + str(readouts[i])]

	all_names = set(names)
	master_names = set(split_names(masters['name'].values))
	new = all_names ^ master_names
	new = list(new)
	print('Number of new flat entries: ',len(new))
	if len(new) > 0:
		new.sort(reverse=True)
		indexer = np.arange(len(new))
		entries = Parallel(n_jobs=num_cores)(delayed(bc_flat_processing)(index,new,flat_list,times,time_frame,save_location,verbose,200,exp_tol) for index in tqdm(indexer, desc='Processing files'))

		entries = pd.concat(entries,ignore_index=True)	
		
		masters = masters.append(entries, ignore_index=True)
		masters.to_csv(package_directory + 'cal_lists/bc_master_flat_list.csv',index=False)

def bc_flat_processing(index,new,flat_list,times,time_frame,save_location,verbose, threshold = 500, exp_tol = 1):
	i = index
	entry = {}
	n = new[i]
	t = int(float(n.split('_')[1]))
	b = n.split('_')[3]
	r = n.split('_')[4]

	tind = (t - times >= 0) & (t - times <= time_frame)
	bind = flat_list['band'] == b
	rind = flat_list['readout'] == r
	ind = tind & bind & rind

	files = flat_list['filename'].values[ind]
	exptimes = flat_list['exptime'].values[ind]
	readouts = flat_list['readout'].values[ind]
	chips = flat_list['chip'].values[ind]
	jds = flat_list['jd'].values[ind]

	master_arr = []
	darks = []
	#try:
	for j in range(len(files)):
		
		hdu = fits.open(files[j])[0]
		header = hdu.header
		data = hdu.data.astype(float)
		data_c = deepcopy(data)
  
		if data_c.shape[1] == 2048:
			for idx in range(4):
				ind = idx*2048//4
				indx = idx*(2048//4 - 1)
				data_c[indx] = np.nan
				data_c[ind] = np.nan
			data_c[-1] = np.nan

		data_collapse = np.nanmean(data_c, axis = 1)
		grads = np.gradient(data_collapse)
		grads[np.isnan(grads)] = 0

		failed_rows = np.where(np.abs(grads) > threshold)[0]

		if len(failed_rows) > 0:
			# print(grads)
			data = data * np.nan

		else:
			saturations = (data > 50000).flatten()
			if sum(saturations) > len(saturations) * 0.3: # if more than 10% of pixels are saturated, set array to nan
				data = data * np.nan
				# print('Saturated')
		
		master_arr.append(data)
		# ZAAACC
		fname, tdiff = get_master_dark(jds[j], exptimes[j], chips[j], readouts[j], 
                                 	   strict=True, exp_tol=exp_tol, telescope = 'bc', shape = data_c.shape[1])
		
		try:
			darks += [fits.open(fname)[0].data]
		except:
			darks += [data * np.nan]
			print('No dark available')

	master_arr = np.array(master_arr)
 
	darks = np.array(darks)
 
	master_arr = master_arr - darks

	mas = np.nanmean(master_arr,axis=0)
	std = np.nanstd(master_arr,axis=0)
	
	header['JDSTART'] = t 
	header['MASTER'] = True
	phdu = fits.PrimaryHDU(data = mas, header = header)
	ehdu = fits.ImageHDU(data = std, header = header)
	hdul = fits.HDUList([phdu, ehdu])

	save_name = save_location + n + '.fits'
	hdul.writeto(save_name,overwrite=True)
	compress = 'gzip -f ' + save_name
	os.system(compress)
	entry['name'] = n
	entry['telescope'] = 'B&C'
	entry['band'] = header['FILTER'].strip()
	entry['exptime'] = header['EXPTIME']
	entry['chip'] = -1
	entry['readout'] = header['READOUTM'].replace(" ", "")
	entry['jd'] = t
	entry['date'] = header['DATE-OBS'].strip()
	entry['filename'] = save_name + '.gz'
	entry['nimages'] = len(master_arr)
	entry['flattype'] = 'sky' # Just assume for now
	entry['shape'] = mas.shape[0]

	if (np.nanmedian(mas) < 15000) | (np.nansum(mas) <= 0):
		print(np.nanmedian(mas), np.nansum(mas))
		note = 'bad'
	else:
		note = 'good'
	entry['note'] = note

	return pd.DataFrame([entry])

# def moa_flat_processing(index,new,flat_list,times,time_frame,save_location,verbose, exp_tol):
# 	i = index
# 	entry = {}
# 	n = new[i]
# 	t = int(n[1:].split('_')[0])
# 	c = int(n.split('_')[-1])
# 	b = n.split('_')[2]

# 	tind = (t - times >= 0) & (t - times <= time_frame)
# 	cind = flat_list['chip'] == c
# 	bind = flat_list['band'] == b
# 	ind = tind & cind & bind

# 	files = flat_list['filename'][ind].tolist()
# 	exptimes = flat_list['exptime'][ind].tolist()

# 	master_arr = []
	# darks = []
	# if len(files) > 10:
	# 	files = files[:10]
	# for j in range(len(files)):
		
	# 	hdu = fits.open(files[j])[0]
	# 	header = hdu.header
	# 	#print(header)
	# 	data = hdu.data.astype(float)

	# 	saturations = (data > 50000).flatten()
	# 	# if more than 10% of pixels are saturated, set array to nan
	# 	if sum(saturations) > len(saturations) * 0.1:
	# 		data = data * np.nan
	# 	master_arr += [data]

	# 	fname, tdiff = get_master_dark(t, exptimes[j], c, readout = -1, exp_tol = exp_tol)
	# 	try:
	# 		darks += [fits.open(fname)[0].data]
	# 	except:
	# 		darks += [data * np.nan]

	# master_arr = np.array(master_arr)
	# darks = np.array(darks)
	# master_arr = master_arr - darks

	# mas = np.nanmean(master_arr,axis=0)
	# std = np.nanstd(master_arr,axis=0)
	# header['JDSTART'] = t 
	# header['MASTER'] = True
	# phdu = fits.PrimaryHDU(data = mas, header = header)
	# ehdu = fits.ImageHDU(data = std, header = header)
	# hdul = fits.HDUList([phdu, ehdu])

	# save_name = save_location + n + '.fits'
	# hdul.writeto(save_name,overwrite=True)
	# compress = 'gzip -f ' + save_name
	
	# with HiddenPrints():
	# 	os.system(compress)
	# entry['name'] = n

	# entry['band'] = header['COLOUR'].strip()
	# entry['chip'] = header['CHIP']
	# entry['exptime'] = header['EXPTIME']
	# entry['jd'] = t
	# entry['date'] = header['DATE-OBS']
	# entry['filename'] = save_name + '.gz'
	# entry['nimages'] = len(master_arr)
	# entry['readout'] = -1
	# entry['telescope'] = 'MOA-II'

	# field = header['FIELD']
	# if 'flat_round' in field:
	# 	flat_type = 'dome'
	# else:
	# 	flat_type = 'sky'
	# entry['field'] = field
	# entry['flat_type'] = flat_type

	# if (np.nanmedian(mas) < 15000) | (np.nansum(mas) <= 0):
	# 	note = 'bad'
	# else:
	# 	if (len(master_arr) < 2) & (flat_type == 'dome'):
	# 		note = 'bad'
	# 	else:
	# 		note = 'good'
	# entry['note'] = note
	
	# field = header['FIELD']
	# if 'flat_round' in field:
	# 	flat_type = 'dome'
	# else:
	# 	flat_type = 'sky'
	# entry['field'] = field
	# entry['flat_type'] = flat_type

def make_master_flats(save_location = '/home/phys/astronomy/rri38/moa/data/master/flat/',redo_bad=False, verbose=False, telescope = 'moa', time_frame=60,num_cores=25, exp_tol = 1):
	if telescope.lower() == 'bc':
		save_location = '/home/phys/astronomy/rri38/fli/data/master/flat/'
		_bc_master_flats(save_location, time_frame, num_cores, verbose, exp_tol)
	else:
		save_location = '/home/phys/astronomy/rri38/moa/data/master/flat/'
		_moa_master_flats(save_location, time_frame, num_cores, verbose, redo_bad, exp_tol)

def _moa_master_flats(save_location, time_frame, num_cores, verbose, redo_bad, exp_tol = 1):
	flat_list = pd.read_csv(package_directory + 'cal_lists/moa_flat_list.csv')
	try:
		masters = pd.read_csv(package_directory + 'cal_lists/moa_master_flat_list.csv')
	except:
		masters = pd.DataFrame(columns=['name','telescope','exptime','jd','date','band','chip','readout','filename','nimages','flattype','note'])
	if redo_bad:
		masters = cut_bad_reductions(masters)
	ind = (flat_list['note'] == 'good')
	flat_list = flat_list[ind]
	flat_list['band'] = flat_list['band'].str.strip()
	times = flat_list['jd'].values.astype(int)
	names = []
	for i in range(len(times)):
		names += [('F' + times.astype(str)[i] + '_' + str(time_frame) + 'd_' + flat_list['band'].values[i]
			 		+ '_' + flat_list['chip'].values.astype(str)[i])]
	all_names = set(names)
	master_names = set(split_names(masters['name'].values))
	new = all_names ^ master_names
	new = list(new)
	print('Number of new flat entries: ',len(new))
	if len(new) > 0:
		new.sort(reverse=True)
		indexer = np.arange(len(new))
		entries = Parallel(n_jobs=num_cores)(delayed(moa_flat_processing)(index,new,flat_list,times,time_frame,save_location,verbose, exp_tol) for index in tqdm(indexer, desc='Processing files'))
		entries = pd.concat(entries,ignore_index=True)	
		
		masters = masters.append(entries, ignore_index=True)
		masters.to_csv(package_directory + 'cal_lists/moa_master_flat_list.csv',index=False)

def moa_flat_processing(index,new,flat_list,times,time_frame,save_location,verbose, exp_tol):
	i = index
	entry = {}
	n = new[i]
	t = int(n[1:].split('_')[0])
	c = int(n.split('_')[-1])
	b = n.split('_')[2]

	tind = (t - times >= 0) & (t - times <= time_frame)
	cind = flat_list['chip'] == c
	bind = flat_list['band'] == b
	ind = tind & cind & bind

	files = flat_list['filename'][ind].tolist()
	exptimes = flat_list['exptime'][ind].tolist()

	master_arr = []
	darks = []
	if len(files) > 10:
		files = files[:10]
	for j in range(len(files)):
		
		hdu = fits.open(files[j])[0]
		header = hdu.header
		#print(header)
		data = hdu.data.astype(float)

		saturations = (data > 50000).flatten()
		# if more than 10% of pixels are saturated, set array to nan
		if sum(saturations) > len(saturations) * 0.1:
			data = data * np.nan
		master_arr += [data]

		fname, tdiff = get_master_dark(t, exptimes[j], c, readout = -1, exp_tol = exp_tol)
		try:
			darks += [fits.open(fname)[0].data]
		except:
			darks += [data * np.nan]

	master_arr = np.array(master_arr)
	darks = np.array(darks)
	master_arr = master_arr - darks

	mas = np.nanmean(master_arr,axis=0)
	std = np.nanstd(master_arr,axis=0)
	header['JDSTART'] = t 
	header['MASTER'] = True
	phdu = fits.PrimaryHDU(data = mas, header = header)
	ehdu = fits.ImageHDU(data = std, header = header)
	hdul = fits.HDUList([phdu, ehdu])

	save_name = save_location + n + '.fits'
	hdul.writeto(save_name,overwrite=True)
	compress = 'gzip -f ' + save_name
	
	with HiddenPrints():
		os.system(compress)
	entry['name'] = n

	entry['band'] = header['COLOUR'].strip()
	entry['chip'] = header['CHIP']
	entry['exptime'] = header['EXPTIME']
	entry['jd'] = t
	entry['date'] = header['DATE-OBS']
	entry['filename'] = save_name + '.gz'
	entry['nimages'] = len(master_arr)
	entry['readout'] = -1
	entry['telescope'] = 'MOA-II'

	field = header['FIELD']
	if 'flat_round' in field:
		flat_type = 'dome'
	else:
		flat_type = 'sky'
	entry['field'] = field
	entry['flat_type'] = flat_type

	if (np.nanmedian(mas) < 15000) | (np.nansum(mas) <= 0):
		print(np.nanmedian(mas), np.nansum(mas))
		note = 'bad'
	else:
		if (len(master_arr) < 2) & (flat_type == 'dome'):
			note = 'bad'
		else:
			note = 'good'
	entry['note'] = note
	
	field = header['FIELD']
	if 'flat_round' in field:
		flat_type = 'dome'
	else:
		flat_type = 'sky'
	entry['field'] = field
	entry['flattype'] = flat_type

	return pd.DataFrame([entry])
	
def make_masters(save_location = '/home/phys/astronomy/rri38/moa/data/master/', time_frame_dark = 4, time_frame_flat = 30, num_cores=25,verbose=False, telescope = 'moa', dark_tolerance = 2):
	make_master_darks(save_location = save_location + 'dark/', time_frame = time_frame_dark, num_cores=num_cores,verbose=verbose, telescope = telescope, dark_tolerance = dark_tolerance)
	# print('!!! Finished darks !!!')
	if verbose:
		print('!!! Finished darks !!!')
	make_master_flats(save_location = save_location + 'flat/', verbose=verbose, telescope = telescope, time_frame=time_frame_flat,num_cores=num_cores, exp_tol = dark_tolerance)
	if verbose:
		print('!!! Finished flats !!!')

if __name__ == '__main__':
	make_masters(verbose=True, num_cores = 1)