import os
from os import path
from astropy.io import fits
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astroquery.astrometry_net import AstrometryNet
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
import astropy.units as u
from astropy.wcs import WCS
from copy import deepcopy
import subprocess
from hidden_prints import HiddenPrints
import gc

from scipy.stats import iqr
from aperture_photom import cal_photom

from scipy.ndimage.filters import convolve
from satellite_detection import sat_streaks

package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'
tmp = os.environ['TMPDIR'] # this is a bash environment variable, add TMPDIR to your .bashrc or equivalent


import warnings
warnings.filterwarnings("ignore")


#class Consume_moa():
	
#	def __init__(self,images,):
		
class pouakai():

	def __init__(self,file,time_tolerence=100,dark_tolerence=10,savepath='',
				 local_astrom=True,verbose=True,rescale=True,plot=False, calibrate=True, center = None, limit_source = None):

		self.verbose = verbose
		self.file = file
		self.savepath = savepath
		self._local_astrom = local_astrom
		self.time_tolerence = time_tolerence
		self.dark_tolerence = dark_tolerence
		self.offset = 500
		# try:
		# 	self.filter_list = np.genfromtxt('/home/phys/astronomy/zgl12/fli_pouakai/pouakai/filters_list.csv', 
		# 					delimiter = ',', dtype = str).flatten().tolist()
		# except:
		self.filter_list = []
		self.fail_flag = ''
		self.rescale = rescale
		self.field_coord = None
		self.plot = plot
		self.center = center
		self.limit_source = limit_source


		self._start_record()
		self._check_dirs()
		self._set_base_name() 
		self._calibrate = calibrate
		
		self._read_science_image()

		self._get_master('dark')
		self._get_master('flat')


		self._setup_fig()
		
		#try:
		self.reduce()
		#except Exception as e:
		#	self.fail_flag = e

		if self.fail_flag != '':
			self._fail_log()
		
		del self
		gc.collect()

	def reduce(self):
		print(self.file)
		#self._check_reduction(reduction)
		self.reduce_image()
		self.save_intermediate()
		
		if self._local_astrom:
			self.wcs_astrometrynet_local()
		else:
			self.wcs_astrometrynet()
		
		self.satellite_search()
		self.Make_mask()
		if self._calibrate:
			self.calculate_zp()
			self.save_fig()
		self.save_intermediate_wcs()
		#self.save_image()
		#print('!!!!!!!!!!!!!!!!!!!!!!!!SAVED')
		#self._update_reduction_log()
		self._record_reduction()
		#print('!!!!!!!!!!!!!!!!!!!!!!!!REDUCTION')
		if self._calibrate:
			#print('!!!!!!!!!!!!!!!!!!!!!!!!Initial Cal')
			self._save_phot_table()
			del self.cal
		del self
	#def _check_reduction(self,reduction):

	def _query_object(self):
		from astroquery.simbad import Simbad
		r = None
		obj = self.header['OBJECT'].strip()
		print(obj)
		try:
			result_table = Simbad.query_object(obj)
			r = result_table.to_pandas()
			c = SkyCoord(r.RA.values[0],r.DEC.values[0], unit=(u.hourangle, u.deg))
			self.field_coord = c
			print('Retrieved with Simbad')
		except:
			print('Simbad failed')
		if r is None:
			try:
				from astroquery.ipac.ned import Ned
				result_table = Ned.query_object(obj)
				r = result_table.to_pandas()
				c = SkyCoord(r.RA.values[0],r.DEC.values[0], unit=(u.deg, u.deg))
				self.field_coord = c
				print('Retrieved with NED')
			except:
				print('NED failed')





	def _fail_log(self):
		document = {'fname': self.file,
					'error':self.fail_flag}
		#error_log = pd.read_csv('cal_lists/error_log.csv')
		self.error_log = document
		#error_log.to_csv('cal_lists/error_log.csv',index=False)

	def _start_record(self):
		"""
		Set up a log to keep track of key steps and parameters.
		"""
		document = {'name':None,'band':None,'telescope':None,
					'exptime':None,'jd':None,'date':None,
					'field':None,'flat':None,
					'dark':None,'tdiff_flat':None,
					'tdiff_dark':None,'zp':None,'zperr':None,
					'maglim5':None,'maglim3':None,'savename':None}
		self.log = document 

	def _update_reduction_log(self):
		new_entry = pd.DataFrame([self.log])

		log = pd.read_csv(package_directory + 'cal_lists/calibrated_image_list.csv')
		log = pd.concat([log, new_entry], ignore_index=True)
		log.to_csv(package_directory + 'cal_lists/calibrated_image_list.csv',index=False)

	def _read_science_image(self):
		"""
		Read in the science image to be calibrated.
		"""
		hdul = fits.open(self.file)
		hdu = hdul[0]
		self.header = hdu.header
		self.raw_image = hdu.data
		self.jd = hdu.header['JD']
		self.filter = hdu.header['FILTER'].strip()
		self.telescope = hdu.header['TELESCOP']
		self.exp_time = int(hdu.header['EXPTIME'])
		try:
			self._field_coords()
		except:
			self._query_object()

		self.log['band'] = self.filter
		self.log['raw_filename'] = self.file
		self.log['jd'] = self.jd
		self.log['telescope'] = self.telescope.strip()
		self.log['exptime'] = self.exp_time
		self.log['date'] = hdu.header['DATE-OBS'].strip()
		self.log['field'] = hdu.header['OBJECT'].strip()
		hdul.close()


	def _field_coords(self):
		"""
		Get the coordinates of the field center. This is not the coordinates for the 
		center of this image.
		"""
		ra = self.header['RA'].strip(' ')
		dec = self.header['DEC'].strip(' ')
		c = SkyCoord(ra,dec, unit=(u.hourangle, u.deg))
		self.field_coord = c
		self.header['RA'] = c.ra.deg
		self.header['DEC'] = c.dec.deg



	def _set_base_name(self):
		"""
		Strip the fluff so that only the base name remains
		"""
		basesname = self.file.split('/')[-1].split('.f')[0].replace(' ','_')
		basename = basesname.split('-')

		bases_list = []

		for i in range(len(basename) - 1):
			bases_list.append(basename[i])
			bases_list.append('_')

		bases_list[-1] = '-' 
		result = ''.join(bases_list)
		filters = basename[-1].split('_')
		self.filter_list.append(filters[-1])

		# Count occurrences of 'g'
		count = sum(1 for s in self.filter_list if s == filters[-1])

		self.base_name = result + str(count).zfill(4) + '_' + filters[-1]

		self.log['name'] = self.base_name
		# np.savetxt('/home/phys/astronomy/zgl12/fli_pouakai/pouakai/filters_list.csv', self.filter_list, delimiter = ',', fmt = '%s')


	def _get_master(self,cal_type):
		"""
		Retrieve the best master calibration image for the science image.
		This cane be used to retrieve flats or darks.
		"""
		if cal_type.lower() == 'flat':
			masters = pd.read_csv('cal_lists/master_flat_list.csv')
			ind = (masters['band'].values == self.filter) & (masters['telescope'].values == self.telescope)
			masters = masters.iloc[ind]

			ind = (masters['note'].values == 'good')
			masters = masters.iloc[ind]

		elif cal_type.lower() == 'dark':
			masters = pd.read_csv('cal_lists/master_dark_list.csv')
			ind = masters['telescope'].values == self.telescope
			masters = masters.iloc[ind]
			exptimes = masters['exptime'].values.astype(int)
			ind = np.where(abs(self.exp_time - exptimes)<=self.dark_tolerence)[0]
			if len(ind) == 0:
				m = 'No master darks with exptime {}'.format(self.exp_time)
				raise ValueError(m)
			masters = masters.iloc[ind]

		else:
			raise ValueError('Only flat and dark are valid options!!')
		file, tdiff = self._find_master(masters)
		if file.split('.')[-1] != 'gz':
			file += '.gz'
		self.log[cal_type] = file
		self.log['tdiff_' + cal_type] = tdiff
		print('tdiff for ' + cal_type + '=' + str(tdiff))
		hdu = fits.open(file)
		data = hdu[0].data
		err =  hdu[1].data
		hdu.close()

		if cal_type.lower() == 'flat':
			self.flat_file = file 
			self.flat = data
			self.flat_err = err
			print('flat',file)
		elif cal_type.lower() == 'dark':
			self.dark_file = file 
			self.dark = data
			self.dark_err = err


	def _find_master(self,masters):
		"""
		Find the best master image for the science image.
		"""
		telescope = self.telescope
		date = self.jd
		tolerence = self.time_tolerence
		telescope_ind = masters['telescope'].values == telescope
		if sum(telescope_ind * 1.) == 0:
			m = f'No master files for telescope {telescope} listed in {masters}'
			raise ValueError(m)
		m = masters.iloc[telescope_ind]
		m_date = m.jd.values

		t_diff = abs(m_date - date)
		t_min = np.nanmin(t_diff)
		if t_min > tolerence:
			m = 'No master file that meets the time tolerence of {}. Closest master is {}'.format(tolerence, t_min)
			raise ValueError(m)
		t_ind = np.argmin(t_diff)

		file = masters['filename'].iloc[t_ind]

		return file, t_diff[t_ind]

	#def _update_header_obj()


	def _update_header_sky(self):
		"""
		Update the fits header to the format needed for photpipe
		"""
		mean, med, std = sigma_clipped_stats(self.image, sigma=3.0)
		self.header['SKYADU'] = (str(np.round(med)), 'median sky')
		self.header['SKYSIG'] = (str(np.round(std)), 'median sky noise')
			
	def _update_header_dark(self):
		self.header['DARKCORR'] = (True,self.dark_file)
			
	def _update_header_flat(self):
		self.header['FLATCORR'] = (True,self.flat_file)
		
	def _update_header_zeropoint(self):
		self.header['ZP'] = (str(np.round(self.zp)),'Calibrimbore ' + self.system)

	def _update_header_satellites(self):
		self.header['SAT'] = (self.sat.satellite,'Satellite in image')
		self.header['SATNUM'] = (self.sat.sat_num,'Number of satellites in image')

	def _update_header_wcs(self):
		#print(self.wcs.wcs.naxis)
		self.header['WCSSOLV'] = 'WCS solved'

		new_header = fits.Header()

		for card in self.wcs.cards:
			if card[0].startswith('CTYPE') or card[0].startswith('CRVAL') or card[0].startswith('CRPIX') or card[0].startswith('CD'):
				self.header[card[0]] = (card[1], card[2], card.comments)

		#self.header.update(new_header.to_header())
		#self.header['CTYPE'] = self.wcs['CTYPE']

	def _check_vars(self):
		"""
		Chech that all images are assigned 
		"""
		m = ''
		if self.flat is None:
			m += 'No flat specified\n'
		if self.dark is None:
			m += 'No dark specified\n'
		
		if self.raw_image is None:
			m += 'No image specified\n'
		if len(m) > 1:
			raise ValueError(m)

	def _check_dirs(self):
		"""
		Check that all reduction directories are constructed
		"""
		dirlist = ['red', 'cal','fig','zp_surface','log','phot_table']
		
		for d in dirlist:
			if not os.path.isdir(self.savepath + d):
				os.mkdir(self.savepath + d)

	def reduce_image(self):
		"""
		Flatten the science image
		"""
		self._check_vars()
		
		image = (self.raw_image - self.dark) / (self.flat/np.nanmedian(self.flat))
		if np.nansum(image) == 0:
			raise ValueError('Image is all NaNs')

		bkg = np.nanmedian(image)
		sub = image - bkg + self.offset

		self.image = sub

		self._add_image(self.raw_image,'A')
		self._add_image(self.flat,'B')
		self._add_image(image,'C')
		self._add_image(image,'D')


	def save_intermediate(self):
		"""
		Save the flattened science image
		"""
		self._update_header_sky()
		self._update_header_dark()
		self._update_header_flat()
		name = self.savepath + 'red/' + self.base_name + '_red.fits'
		name = name.replace(' ','_')
		self.red_name = name + '.gz'


		if self.verbose:
			print('Saving intermediated calibrated file')
		fits.writeto(name,self.image,header=self.header,overwrite=True)
		compress = 'gzip -f ' + name
		os.system(compress)


	def save_image(self):
		"""
		Save the final calibrated science image
		"""
		self._update_header_sky()
		self._update_header_dark()
		self._update_header_flat()
		self._update_header_wcs
		name = self.savepath + 'cal/' + self.base_name + '_cal.fits'
		name = name.replace(' ','_')
		self.cal_name = name + '.gz'
		phdu = fits.PrimaryHDU(data = self.image, header = self.header)
		mhdu = fits.ImageHDU(data = self.mask, header = self.header)
		hdul = fits.HDUList([phdu, mhdu])
		print(hdul)
		if self.verbose:
			print('Saving final calibrated image')
		
		hdul.writeto(name,overwrite=True)
		compress = 'gzip -f ' + name
		os.system(compress)
		self.log['savename'] = self.cal_name

	def wcs_astrometrynet(self,timeout=200):
		"""
		Calculate the image wcs using the portal for astrometry.net
		"""
		ast = AstrometryNet()
		ast.api_key = 'csffdfuichpbiata'
		attempt = 0
		solved = False

		attempts_allowed = 16

		#while (attempt < 10) & (not solved):

		# print(self.center, type(self.center))


		for attempt in range(attempts_allowed):

			try:
				with HiddenPrints():
					if (isinstance(self.center, list)) or (isinstance(self.center, tuple)) or (isinstance(self.center, np.ndarray)):
						print('Using given coordinates')
						wcs_head = ast.solve_from_image(self.file, solve_timeout = timeout, 
													scale_units = 'arcminwidth', scale_type='ul', 
													scale_upper = 26, scale_lower = 24, 
													center_ra = self.center[0], center_dec = self.center[1], 
													radius = 0.5, crpix_center=True)
					
					else:
						print('fail_print')
						wcs_head = ast.solve_from_image(self.file, solve_timeout = timeout, 
													scale_units = 'arcminwidth', scale_type='ul', 
													scale_upper = 26, scale_lower = 24, crpix_center=True)

					# else:
					# 	wcs_head = ast.solve_from_image(self.file, solve_timeout = timeout, 
					# 									scale_units = 'arcminwidth', scale_type='ul', 
					# 									scale_upper = 26, scale_lower = 24)

					del wcs_head['COMMENT']
					del wcs_head['HISTORY']
					solved = True
			except:
				wcs_head = 0

			if wcs_head != 0:
				solved = True
				break

		self._wcs_solution = True		
		if solved == False:
			self._wcs_solution = False
			raise ValueError('Could not solve WCS in {} attempts'.format(attempt + 1))

		new_head = deepcopy(self.header)
		for key in wcs_head:
			new_head[key] = (wcs_head[key],wcs_head.comments[key])

		if self.verbose:
			print('Solved WCS, saving file')
		self.header = new_head
		self.wcs = WCS(self.header)

	def wcs_astrometrynet_local(self):
		"""
		Calculate the image wcs using the local libraries for astrometry.net
		"""
		# a reasonable search radius is already selected (2deg)
		save_path = 'wcs_tmp/' + self.base_name + '/'
		real_save_path = self.savepath + 'red/' + save_path
		#if not path.exits(real_save_path):
		os.mkdir(real_save_path)
	
		name = save_path + self.base_name + '_wcs'
		real_name = real_save_path + self.base_name + '_wcs'

		if self.field_coord is not None:
			astrom_call = f"solve-field --no-plots --scale-units arcminwidth --scale-low 24 --scale-high 26 --temp-dir {tmp} -O -o {name} -p --ra {self.field_coord.ra.deg} --dec {self.field_coord.dec.deg} --radius 0.4 {self.red_name}"
		else:
			astrom_call = f"solve-field --no-plots --scale-units arcminwidth --scale-low 24 --scale-high 26 --temp-dir {tmp} -O -o {name} -p {self.red_name}"
		try:
			os.system(astrom_call)
			#try:
			wcs_header = fits.open(real_name + '.new')[0].header
			# get rid of all the astrometry.net junk in the header 
			del wcs_header['COMMENT']
			del wcs_header['HISTORY']
			self.header = wcs_header
			self.wcs = WCS(self.header)

			if self.verbose:
				print('Solved WCS')
			
			clear = 'rm -rf ' + real_save_path
			os.system(clear)

			self._wcs_solution = True
		except:
			self._wcs_solution = False

		if self.verbose:
			print('WCS tmp files cleared')


	def save_intermediate_wcs(self):
		"""
		Save the intermediate image with a wcs solution.
		"""
		if not os.path.exists(self.savepath + 'cal/'):
			os.makedirs(self.savepath + 'cal/')

		name = self.savepath + 'cal/' + self.base_name + '_wcs.fits'
		name = name.replace(' ','_')
		self.cal_name = name + '.gz'
		self.wcs_name = name

		if self.verbose:
			print('Saving intermediated wcs file')
		fits.writeto(name,self.image,header=self.header,overwrite=True)

		compress = 'gzip -f ' + name
		os.system(compress)
		self.log['savename'] = self.cal_name


	def calculate_zp(self,threshold=3,model='ckmodel'):
		"""
		Use calibrimbore to calculate the zeropoint for the image and the magnitude limits.
		"""
		if self.verbose:
			print('Calculating zeropoint')

		if self.log['exptime'] < (2.5 * 60):
			brightlim = 5
		else:
			brightlim = 10
		mask = ((self.mask & 2) + (self.mask & 4) + (self.mask & 8) + (self.mask & 16))
		mask[mask > 0] = 1
		if self.plot:
			self.cal = cal_photom(data=self.image,wcs=self.wcs,mask=mask, header=self.header,
								threshold=threshold,cal_model=model,ax=self.fig_axis['F'],
								brightlim=brightlim,rescale=self.rescale,limit_source=self.limit_source)
		else:
			self.cal = cal_photom(data=self.image,wcs=self.wcs,mask=mask, header=self.header,
								threshold=threshold,cal_model=model,plot=False,
								brightlim=brightlim,rescale=self.rescale,limit_source=self.limit_source)
      
		#self._add_image(self.cal.zp_surface,'E',colorbar=True)
		self._add_image(self.cal.data,'C')
		self._add_satellite_trail('C')
		self.image = self.cal.data

		self.header['ZP'] = (str(np.round(self.cal.zp,4)), 'Calibrimbore zeropoint')
		self.header['ZPERR'] = (str(np.round(self.cal.zp_std,4)), 'Calibrimbore zeropoint error')
		self.header['MAGLIM5'] = (str(np.round(self.cal.maglim5)), '5 sig mag lim')
		self.header['MAGLIM3'] = (str(np.round(self.cal.maglim3)), '3 sig mag lim')
		self.log['zp'] = np.round(self.cal.zp,4)
		self.log['zperr'] = np.round(self.cal.zp_std,4)
		self.log['maglim5'] = np.round(self.cal.maglim5,4)
		self.log['maglim3'] = np.round(self.cal.maglim3,4)

		ind = np.isfinite(self.cal.zps)
		if self.plot:
			self.fig_axis['D'].plot(self.cal.source_x[ind],self.cal.source_y[ind],'r+')
			self._zp_hist()
			#self._zp_color()
			self.cal.mag_limit_fig(self.fig_axis['F'])
		self._save_zp_surface()
		if self.verbose:
			print('Zeropoint found to be ' + str(np.round(self.cal.zp,2)))

	def _save_zp_surface(self):
		path = f'{self.savepath}/zp_surface/{self.base_name}_zp_surface'
		np.save(path,self.cal.zp_surface)

	def _zp_hist(self):
		"""
		Create a histogram of the zeropoint distribution
		"""
		if self.plot:
			zps = self.cal.zps#[self.cal.good]
			zps = zps[np.isfinite(zps)]
			#b = int(abs(np.nanmax(zps) - np.nanmin(zps) /(2*iqr(zps)*len(zps)**(-1/3))))
			self.fig_axis['E'].hist(zps,alpha=0.5)		
			med = self.cal.zp
			high = self.cal.zp+self.cal.zp_std
			low = self.cal.zp-self.cal.zp_std
			self.fig_axis['E'].axvline(med,color='k',ls='--')
			self.fig_axis['E'].axvline(low,color='k',ls=':')
			self.fig_axis['E'].axvline(high,color='k',ls=':')
			s = ('$zp='+str((np.round(med,2)))+'^{+' + 
				str(np.round(high-med,2))+'}_{'+
				str(np.round(low-med,2))+'}$')
			self.fig_axis['E'].annotate(s,(0.7,0.8),fontsize=10,xycoords='axes fraction')
			self.fig_axis['E'].set_xlabel('zeropoint',fontsize=15)
			self.fig_axis['E'].set_ylabel('Occurrence',fontsize=15)

	def _zp_color(self):
		"""
		Create a figure showing the zeropoint evolving with color.
		"""
		zps = self.cal.zps
		ind = np.isfinite(zps)
		gr = (self.cal.sauron.cat_mags['g'] - self.cal.sauron.cat_mags['r']).values[ind]

		self.fig_axis['H'].plot(gr,zps[ind],'.')
		self.fig_axis['H'].set_ylabel('zeropoint',fontsize=15)
		self.fig_axis['H'].set_xlabel('$g-r$',fontsize=15)	
		med = self.cal.zp
		high = self.cal.zp+self.cal.zp_std
		low = self.cal.zp-self.cal.zp_std
		self.fig_axis['H'].axhline(med,color='k',ls='--')
		self.fig_axis['H'].axhline(low,color='k',ls=':')
		self.fig_axis['H'].axhline(high,color='k',ls=':')


	def save_fig(self):
		"""
		Save a diagnostic figure
		"""
		if self.plot:
			name = self.savepath + 'fig/' + self.base_name + '_diag.png'
			#self.fig.set_tight_layout(True)
			self.fig.savefig(name)

	def _setup_fig(self):
		"""
		Set up the large diagnostic plot figure
		"""
		if self.plot:
			self.fig = plt.figure(figsize=(8.27,11.69),constrained_layout=True)
			self.fig_axis = self.fig.subplot_mosaic(
												"""
												ABC
												ABC
												DEF
												"""
											)
			self.fig_axis['A'].set_title('Raw image',fontsize=15)
			self.fig_axis['B'].set_title('Flat image',fontsize=15)
			self.fig_axis['C'].set_title('Reduced image',fontsize=15)
			self.fig_axis['D'].set_title('Calibration sources',fontsize=15)
			self.fig_axis['E'].set_title('Zeropoint distribution',fontsize=15)
			#self.fig_axis['F'].set_title('Zeropoint colour',fontsize=15)
			self.fig_axis['F'].set_title('Signal-Noise Limit',fontsize=15)

	def _add_image(self,image,ax_ind,colorbar=False):
		"""
		Add the provided image to the provided axis.
		"""
		if self.plot:
			vmin = np.percentile(image,16)
			vmax = np.percentile(image,84)
			im = self.fig_axis[ax_ind].imshow(image,origin='lower',
										vmin=vmin,vmax=vmax)
			if colorbar:
				self.fig.colorbar(im,ax=self.fig_axis[ax_ind],fraction=0.046, pad=0.04)
	
	def _add_satellite_trail(self,ax_ind):
		if self.plot:
			for consolidated_line in self.sat.consolidated_lines:
				angle = consolidated_line[0]
				points = np.array(consolidated_line[1])
				x = points[:,0]
				y = points[:,1]
				coefs = np.polyfit(x, y, 1)
				slope = coefs[0]
				intercept = coefs[1]
				x1 = np.min(x)
				y1 = int(slope * x1 + intercept)
				x2 = np.max(x)
				y2 = int(slope * x2 + intercept)
				# create the x and y points for the consolidated line
				x_points = [x1, x2]
				y_points = [y1, y2]

				# plot the consolidated line
				self.fig_axis[ax_ind].plot(x_points, y_points, '--r', label='Satellite',alpha = .5)
				self.fig_axis[ax_ind].legend(loc='upper left')


	def _record_reduction(self):
		"""
		Save the log to the global reduction log file.
		"""
		#log = pd.read_csv('cal_lists/calibrated_image_list.csv')
		new_entry = pd.DataFrame([self.log])
		#log = pd.concat([log, new_entry], ignore_index=True)
		new_entry.to_csv(f'{self.savepath}log/{self.base_name}.csv')
		#log.to_csv('cal_lists/calibrated_image_list.csv',index=False)



	def _flat_mask(self,lowlim=0.8,highlim=1.2,buffer=3):
		mask = deepcopy(self.flat) / np.nanmedian(deepcopy(self.flat))
		mask = (mask < lowlim) | (mask > highlim)
		kernel = np.ones((buffer,buffer))
		mask = convolve(mask,kernel)
		mask = mask.astype(int)
		return mask

	def _saturaton_mask(self,satlimit=3.5e4,buffer=3):
		"""
		An agressive limit is set here to catch overflow pixels.
		"""
		mask = deepcopy(self.raw_image)
		mask = mask > satlimit
		mask = convolve(mask,np.ones((buffer,buffer)))
		mask = mask.astype(int)
		return mask

	def _load_bad_pix_mask(self):
		bpix = np.load(f'badpix/telescope_{self.telescope}_bpix.npy')
		return bpix.astype(int)

	def Make_mask(self):
		
		saturation_mask = self._saturaton_mask() * 2
		flat_mask = self._flat_mask() * 4
		#satellite_mask = self.sat.total_mask * 8
		#bpix = self._load_bad_pix_mask() * 16

		self.mask = flat_mask | saturation_mask #| satellite_mask #| bpix

		self._update_header_mask_bits()


	def _update_header_mask_bits(self):
		'''
		Function to update the header mask bits.
		'''
		#head['STARBIT']  = (1, 'bit value for normal sources')
		self.header['SATBIT']   = (2, 'bit value for saturated sources')
		self.header['FLATBIT'] = (4, 'bit value for bad flat')
		#self.header['SATBIT'] = (8, 'bit value for satellite pixels')
		self.header['USERBIT']  = (16, 'Bad pixels')
		#head['SNBIT']    = (32, 'bit value for SN list')


	def satellite_search(self):
		self.sat = sat_streaks(self.image,run=True)
		self._update_header_satellites()
		

	def _save_phot_table(self):
		name = self.savepath + 'phot_table/' + self.base_name + '_phot.fits'
		ind = self.cal.ap_photom.mag.values < self.cal.maglim3
		tab = self.cal.ap_photom.iloc[ind]
		ra,dec = self.wcs.all_pix2world(tab['xcenter'].values,tab['ycenter'].values,0)
		tab['ra'] = ra
		tab['dec'] = dec
		rec = np.rec.array([tab['gaiaID'].values,tab['xcenter'].values,tab['ycenter'].values,tab['ra'].values,tab['dec'].values,tab['counts'].values,tab['e_counts'].values,tab['mag'].values,tab['e_mag'].values,tab['snr'].values],
							formats='int,float32,float32,float32,float32,float32,float32,float32,float32,float32',
							names='gaiaID,xcenter,ycenter,ra,dec,counts,counts_e,mag,mag_e,snr' )

		hdu = fits.BinTableHDU(data=rec,header=self.header)
		
		hdu.writeto(name,overwrite=True)
