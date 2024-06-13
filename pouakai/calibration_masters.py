from astropy.io import fits
import pandas as pd
import numpy as np
from glob import glob
import os
from copy import deepcopy
from joblib import Parallel, delayed

def split_names(files):
	names = [x.split('-')[0] for x in files]
	return names


def make_master_darks(save_location = '/home/phys/astronomy/rri38/fli/data/master/dark/',time_frame=1,num_cores=25,verbose=False):
	# make save_location an environment variable
	dark_list = pd.read_csv('cal_lists/dark_list.csv')
	masters = pd.read_csv('cal_lists/master_dark_list.csv')
	master_name = assign_master_name(dark_list)
	all_names = set(master_name)
	master_names = set(masters['name'].values)

	new = all_names ^ master_names
	new = list(new)
	print('Number of new dark entries: ',len(new))
	if len(new) > 0:
		new.sort(reverse=True)
		indexer = np.arange(len(new))
		entries = Parallel(n_jobs=num_cores)(delayed(dark_processing)(index,new,dark_list,time_frame,save_location,verbose)  for index in indexer)
	for entry in entries:
		masters = pd.concat([masters, entry],ignore_index=True)
	masters.to_csv('cal_lists/master_dark_list.csv',index=False)

def assign_master_name(darks):
	jd = darks['jd'].values.astype(int).astype(str)
	exptime = darks['exptime'].values.astype(int).astype(str)
	master_name = []
	for i in range(len(jd)):
		master_name += [f'dark_{jd[i]}_{exptime[i]}s']
	return master_name

def dark_processing(index,names,dark_list,time_frame,save_location,verbose):
	entry = {}
	name = names[index]
	t = int(name.split('_')[1])
	exptime = int(name.split('_')[2].split('s')[0])
	times = dark_list['jd'].values.astype(int)
	tind = (t - times >= -time_frame/2) & (t - times <= time_frame)
	expind = dark_list['exptime'].values.astype(int) == exptime
	ind = tind & expind

	files = dark_list['filename'].values[ind]
	master = []
	try:
		for file in files:
			hdu = fits.open(file)[0]
			header = hdu.header
			data = hdu.data
			master += [data]
		master = np.array(master)
		#print('made array')
		if verbose:
			print('Used ',len(master),' images in median')
		m = np.nanmedian(master,axis=0)
		std = np.nanstd(master,axis=0)
		time = np.nanmean(dark_list['jd'].iloc[ind])
		#print('calc mean')
		header['JDSTART'] = time 
		header['MASTER'] = True
		phdu = fits.PrimaryHDU(data = m, header = header)
		ehdu = fits.ImageHDU(data = std, header = header)
		hdul = fits.HDUList([phdu, ehdu])

		save_name = save_location + name + '.fits'
		print('saving')
		hdul.writeto(save_name,overwrite=True)

		compress = 'gzip -f ' + save_name
		os.system(compress)
		print('saved')
		entry['name'] = name
		entry['telescope'] = dark_list['telescope'].values[ind][0]
		entry['exptime'] = header['EXPTIME']
		entry['jd'] = time
		entry['date'] = header['DATE-OBS']
		entry['nimages'] = len(master)
		entry['filename'] = save_name + '.gz'
		if len(master) < 3:
			note = 'bad'
		else:
			note = 'good'
		entry['note'] = note
		if verbose:
			print('Done ', name)
		entry = pd.DataFrame([entry])
			
		return entry
	except Exception as e:
		print(e)
		print('Something went wrong...')

def get_master_dark(jd,exptime,strict=True,tol=10,exp_tol=5):
	"""
	ytdhgvj
	"""
	darks = pd.read_csv('cal_lists/master_dark_list.csv')
	if strict:
		ind = darks['note'].values == 'good'
		darks = darks.iloc[ind]
	
	dexptime = darks['exptime'].values
	exp_ind = abs(dexptime.astype(int) - int(exptime)) < exp_tol
	good = darks.iloc[exp_ind]

	if len(good) > 0:
		djd = good['jd'].values
		diff = jd - djd
		min_ind = np.argmin(abs(diff))
		t_diff = diff[min_ind]
		dark = good.iloc[min_ind]
		fname = dark['filename']
		#print('flat exp:{}, chip:{}'.format(exptime,chip))
		if abs(t_diff) < tol:
			return fname, t_diff
		else:
			return 'none', -999	
	else:
		print('no exposure')
		print(int(exptime))
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


def make_master_flats(save_location = '/home/phys/astronomy/rri38/fli/data/master/flat/',time_frame=60,num_cores=25, verbose=False):
	# make save_location an environment variable
	flat_list = pd.read_csv('cal_lists/flat_list.csv')
	masters = pd.read_csv('cal_lists/master_flat_list.csv')
	ind = (flat_list['note'].values == 'good')
	flat_list = flat_list.iloc[ind]
	flat_list['band'] = flat_list['band'].str.strip()
	times = flat_list['jd'].values.astype(int)
	names = []
	for i in range(len(times)):
		names += [('flat_' + times.astype(str)[i] + '_' + str(time_frame) + 'd_' + flat_list['band'].values[i])]

	all_names = set(names)
	master_names = set(split_names(masters['name'].values))
	new = all_names ^ master_names
	new = list(new)
	print('Number of new flat entries: ',len(new))
	if len(new) > 0:
		new.sort(reverse=True)
		indexer = np.arange(len(new))
		entries = Parallel(n_jobs=num_cores)(delayed(flat_processing)(index,new,flat_list,times,time_frame,save_location,verbose) for index in indexer)
		entries = pd.concat(entries,ignore_index=True)	
		
		masters = pd.concat([masters, entries], ignore_index=True)
		masters.to_csv('cal_lists/master_flat_list.csv',index=False)

def flat_processing(index,new,flat_list,times,time_frame,save_location,verbose, threshold = 350):
	i = index
	entry = {}
	n = new[i]
	t = int(n.split('_')[1])
	b = n.split('_')[3]

	tind = (t - times >= 0) & (t - times <= time_frame)
	bind = flat_list['band'].values == b
	ind = tind & bind

	files = flat_list['filename'].values[ind]
	exptimes = flat_list['exptime'].values[ind]
	if verbose:
		print('Num flats: ', len(files))

	master_arr = []
	darks = []
	#try:
	for j in range(len(files)):
		
		hdu = fits.open(files[j])[0]
		header = hdu.header
		#print(header)
		data = hdu.data.astype(float)
		data_c = deepcopy(data)
		for idx in range(4):
			ind = idx*512
			indx = idx*511
			data_c[indx] = np.nan
			data_c[ind] = np.nan
		data_c[-1] = np.nan

		data_collapse = np.nanmean(data_c, axis = 1)
		grads = np.gradient(data_collapse)

		failed_rows = np.where(abs(grads) > threshold)[0]

		if len(failed_rows) > 0:
			print('image ', files[j], 'has readout error')
			data = data * np.nan
		
		else:
			saturations = (data > 50000).flatten()
			# if more than 10% of pixels are saturated, set array to nan
			if sum(saturations) > len(saturations) * 0.1:
				print('image ', files[j], ' is saturated')
				data = data * np.nan
		
		master_arr += [data]

		fname, tdiff = get_master_dark(t, exptimes[j])
		
		try:
			darks += [fits.open(fname)[0].data]
		except:
			darks += [data * 0]
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
	print(f'saving: {save_name}')
	hdul.writeto(save_name,overwrite=True)
	compress = 'gzip -f ' + save_name
	os.system(compress)
	print('saved')
	entry['name'] = n
	entry['telescope'] = header['TELESCOP'].strip()
	entry['band'] = header['FILTER'].strip()
	entry['exptime'] = header['EXPTIME']
	entry['jd'] = t
	entry['date'] = header['DATE-OBS'].strip()
	entry['filename'] = save_name + '.gz'
	entry['nimages'] = len(master_arr)

	if (np.nanmedian(mas) < 15000) | (np.nansum(mas) <= 0):
		note = 'bad'
	else:
		note = 'good'
	entry['note'] = note


	if verbose:
		print('Done ', n)
	return pd.DataFrame([entry])
	#except:
	#	print('something went wrong...')

def make_masters(verbose=True):
	make_master_darks(verbose=True)
	if verbose:
		print('!!! Finished darks !!!')
	make_master_flats(verbose=True)

	if verbose:
		print('!!! Finished flats !!!')


if __name__ == '__main__':
	make_master_darks(verbose=True)
	print('!!! Finished darks !!!')
	make_master_flats(verbose=True)

	print('!!! Finished flats !!!')