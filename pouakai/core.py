import os
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

from scipy.stats import iqr
from aperture_photom import ap_photom

from scipy.ndimage.filters import convolve

import warnings
warnings.filterwarnings("ignore")


#class Consume_moa():
	
#	def __init__(self,images,):
		
class pouakai():

	def __init__(self,file,reduction='full',time_tolerence=100,dark_tolerence=10,savepath='',
				 local_astrom=True,verbose=True,update_cals=False):

		#self._update_cals(update_cals)
		self.verbose = verbose
		self.file = file 
		self.savepath = savepath
		self._local_astrom = local_astrom
		self.time_tolerence = time_tolerence
		self.dark_tolerence = dark_tolerence
		self.offset = 500

		self._start_record()
		self._check_dirs()
		self._set_base_name() 
		
		self._read_science_image()

		self._get_master('dark')
		self._get_master('flat')


		self._setup_fig()
		

		self.reduce(reduction)


	def reduce(self,reduction):
		#self._check_reduction(reduction)
		self.reduce_image()
		self.save_intermediate()
		if self._local_astrom:
			self.wcs_astrometrynet_local()
		else:
			self.wcs_astrometrynet()

		self.Make_mask()
		self.calculate_zp()
		self.save_fig()
		self.save_image()
		self._record_reduction()
	#def _check_reduction(self,reduction):



	def _start_record(self):
		"""
		Set up a log to keep track of key steps and parameters.
		"""
		document = {'name':None,'band':None,'chip':None,
					'exptime':None,'jd':None,'date':None,
					'field':None,'filename':None,'flat':None,
					'dark':None,'tdiff_flat':None,
					'tdiff_dark':None,'zp':None,'zperr':None,
					'maglim5':None,'maglim3':None,'savename':None}
		self.log = document 


	def _read_science_image(self):
		"""
		Read in the science image to be calibrated.
		"""
		hdu = fits.open(self.file)[0]
		self.header = hdu.header
		self.raw_image = hdu.data
		self.jd = hdu.header['JDSTART']
		self.filter = hdu.header['COLOUR']
		self.chip = hdu.header['CHIP']
		self.exp_time = hdu.header['EXPTIME']
		self._field_coords()

		self.log['band'] = self.filter
		self.log['raw_filename'] = self.file
		self.log['jd'] = self.jd
		self.log['chip'] = self.chip
		self.log['exptime'] = self.exp_time
		self.log['date'] = hdu.header['DATE-OBS']
		self.log['field'] = hdu.header['FIELD']


	def _field_coords(self):
		"""
		Get the coordinates of the field center. This is not the coordinates for the 
		center of this image.
		"""
		ra = self.header['RA'].strip(' ')
		dec = self.header['DEC'].strip(' ')
		c = SkyCoord(ra,dec, unit=(u.hourangle, u.deg))
		self.field_coord = c



	def _set_base_name(self):
		"""
		Strip the fluff so that only the base name remains
		"""
		self.base_name = self.file.split('/')[-1].split('.f')[0]
		self.log['name'] = self.base_name

	def _get_master(self,cal_type):
		"""
		Retrieve the best master calibration image for the science image.
		This cane be used to retrieve flats or darks.
		"""
		if cal_type.lower() == 'flat':
			masters = pd.read_csv('cal_lists/master_flat_list.csv')
			ind = (masters['band'].values == self.filter) & (masters['chip'].values == self.chip)
			masters = masters.iloc[ind]

			ind = (masters['note'].values == 'good') & (masters['flat_type'].values == 'dome')
			masters = masters.iloc[ind]

		elif cal_type.lower() == 'dark':
			masters = pd.read_csv('cal_lists/master_dark_list.csv')
			ind = masters['chip'].values == self.chip
			masters = masters.iloc[ind]
			exptimes = masters['exptime'].values
			ind = np.where(abs(self.exp_time - exptimes)<self.dark_tolerence)[0]
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
		chip = self.chip
		date = self.jd
		tolerence = self.time_tolerence

		chip_ind = masters['chip'].values == chip
		if len(chip_ind) == 0:
			m = 'No master files for chip {} listed in {}'.format(chip,masters)
			raise ValueError(m)
		m = masters.iloc[chip_ind]
		m_date = m.jd.values

		t_diff = abs(m_date - date)
		t_min = np.nanmin(t_diff)
		if t_min > tolerence:
			m = 'No master file in {} that meets the time tolerence of {}'.format(masters, tolerence)
			raise ValueError(m)
		t_ind = np.argmin(t_diff)

		file = masters['filename'].iloc[t_ind]

		return file, t_diff[t_ind]




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
		dirlist = ['wcs','red','red/wcs_tmp','cal','fig']
		for d in dirlist:
			if not os.path.isdir(self.savepath + d):
				os.mkdir(self.savepath + d)

	def reduce_image(self):
		"""
		Flatten the science image
		"""
		self._check_vars()
		
		image = (self.raw_image - self.dark) / (self.flat/np.nanmedian(self.flat))
		print('raw ', np.nanmean(self.raw_image))
		print('dark ', np.nanmean(self.dark))
		print('flat ', np.nanmean(self.flat))
		print('image ', np.nanmean(image))
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
		name = self.savepath + 'cal/' + self.base_name + '_cal.fits'
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

	def wcs_astrometrynet(self,timeout=120):
		"""
		Calculate the image wcs using the portal for astrometry.net
		"""
		ast = AstrometryNet()
		ast.api_key = 'csffdfuichpbiata'
		attempt = 0
		solved = False
		while (attempt < 10) & (not solved):
			try:
				wcs_head = ast.solve_from_image(self.file,solve_timeout=timeout)
				del wcs_head['COMMENT']
				del wcs_head['HISTORY']
				solved = True
			except:
				attempt += 1
		if (attempt > 10) | (not solved):
			raise ValueError('Could not solve WCS in {} attempts'.format(attempt))
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
		astrom_call = "solve-field -O -o {savename} -p --ra {ra} --dec {dec} --radius 2 {file}"

		save_path = 'wcs_tmp/' + self.base_name + '/'
		real_save_path = self.savepath + 'red/' + save_path
		os.mkdir(real_save_path)
	
		name = save_path + self.base_name + '_wcs'
		real_name = real_save_path + self.base_name + '_wcs'
		print('!!!',name)
		print(save_path)
		print(self.base_name)
		solver = astrom_call.format(savename = name, ra = self.field_coord.ra.deg,
									dec = self.field_coord.dec.deg, file = self.red_name)
		os.system(solver)

		wcs_header = fits.open(real_name + '.new')[0].header
		# get rid of all the astrometry.net junk in the header 
		del wcs_header['COMMENT']
		del wcs_header['HISTORY']
		self.header = wcs_header
		self.wcs = WCS(self.header)

		if self.verbose:
			print('Solved WCS')
		
		clear = 'rm -r ' + real_save_path
		#os.system(clear)

		if self.verbose:
			print('WCS tmp files cleared')


	def save_intermediate_wcs(self):
		"""
		Save the intermediate image with a wcs solution.
		"""
		name = self.savepath + 'wcs/' + self.base_name + '_wcs.fits'
		self.wcs_name = name

		if self.verbose:
			print('Saving intermediated wcs file')
		fits.writeto(name,self.image,header=self.header,overwrite=True)

	def calculate_zp(self,threshold=10,model='ckmodel'):
		"""
		Use calibrimbore to calculate the zeropoint for the image and the magnitude limits.
		"""
		if self.verbose:
			print('Calculating zeropoint')

		if self.log['exptime'] < (2.5 * 60):
			brightlim = 13
		else:
			brightlim = 15

		self.cal = ap_photom(data=self.image,wcs=self.wcs,mask=self.mask, header=self.header,
							 threshold=threshold,cal_model=model,ax=self.fig_axis['I'],
							 brightlim=brightlim)
		self._add_image(self.cal.zp_surface,'E',colorbar=True)
		self._add_image(self.cal.data,'F')
		self.image = self.cal.data
		self.header['ZP'] = (str(np.round(self.cal.zp,2)), 'Calibrimbore zeropoint')
		self.header['ZPERR'] = (str(np.round(self.cal.zp_std,2)), 'Calibrimbore zeropoint error')
		self.header['MAGLIM5'] = (str(np.round(self.cal.maglim5)), '5 sig mag lim')
		self.header['MAGLIM3'] = (str(np.round(self.cal.maglim3)), '3 sig mag lim')
		self.log['zp'] = self.cal.zp
		self.log['zperr'] = self.cal.zp_std
		self.log['maglim5'] = self.cal.maglim5
		self.log['maglim3'] = self.cal.maglim3

		self.fig_axis['D'].plot(self.cal.source_x[self.cal.good],self.cal.source_y[self.cal.good],'r.')
		self._zp_hist()
		#self._zp_color()
		self.cal.mag_limit_fig(self.fig_axis['I'])
		if self.verbose:
			print('Zeropoint found to be ' + str(np.round(self.cal.zp,2)))

	def _zp_hist(self):
		"""
		Create a histogram of the zeropoint distribution
		"""
		zps = self.cal.zps
		zps = zps[np.isfinite(zps)]
		#b = int(abs(np.nanmax(zps) - np.nanmin(zps) /(2*iqr(zps)*len(zps)**(-1/3))))
		self.fig_axis['G'].hist(zps,alpha=0.5)		
		med = self.cal.zp
		high = self.cal.zp+self.cal.zp_std
		low = self.cal.zp-self.cal.zp_std
		self.fig_axis['G'].axvline(med,color='k',ls='--')
		self.fig_axis['G'].axvline(low,color='k',ls=':')
		self.fig_axis['G'].axvline(high,color='k',ls=':')
		s = ('$zp='+str((np.round(med,2)))+'^{+' + 
			str(np.round(high-med,2))+'}_{'+
			str(np.round(low-med,2))+'}$')
		self.fig_axis['G'].annotate(s,(.7,.8),fontsize=10,xycoords='axes fraction')
		self.fig_axis['G'].set_xlabel('zeropoint',fontsize=15)
		self.fig_axis['G'].set_ylabel('Occurrence',fontsize=15)

	def _zp_color(self):
		"""
		Create a figure showing the zeropoint evolving with color.
		"""
		zps = self.cal.zps[self.cal.good]
		gr = (self.cal.cat_mags['g'] - self.cal.cat_mags['r']).values[self.cal.good]

		self.fig_axis['H'].plot(gr,zps,'.')
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
		name = self.savepath + 'fig/' + self.base_name + '_diag.pdf'
		#self.fig.set_tight_layout(True)
		self.fig.savefig(name)

	def _setup_fig(self):
		"""
		Set up the large diagnostic plot figure
		"""
		self.fig = plt.figure(figsize=(8.27,11.69),constrained_layout=True)
		self.fig_axis = self.fig.subplot_mosaic(
											"""
											ABC
											ABC
											DEF
											DEF
											GHI
											"""
										  )
		self.fig_axis['A'].set_title('Raw image',fontsize=15)
		self.fig_axis['B'].set_title('Flat image',fontsize=15)
		self.fig_axis['C'].set_title('Reduced image',fontsize=15)
		self.fig_axis['D'].set_title('Calibration sources',fontsize=15)
		self.fig_axis['E'].set_title('Zeropoint correction',fontsize=15)
		self.fig_axis['F'].set_title('Rescaled image',fontsize=15)
		self.fig_axis['G'].set_title('Zeropoint distribution',fontsize=15)
		self.fig_axis['H'].set_title('Zeropoint colour',fontsize=15)
		self.fig_axis['I'].set_title('Signal-Noise Limit',fontsize=15)

	def _add_image(self,image,ax_ind,colorbar=False):
		"""
		Add the provided image to the provided axis.
		"""
		vmin = np.percentile(image,16)
		vmax = np.percentile(image,84)
		im = self.fig_axis[ax_ind].imshow(image,origin='lower',
									 vmin=vmin,vmax=vmax)
		if colorbar:
			self.fig.colorbar(im,ax=self.fig_axis[ax_ind],fraction=0.046, pad=0.04)
			

	def _record_reduction(self):
		"""
		Save the log to the global reduction log file.
		"""
		log = pd.read_csv('cal_lists/calibrated_image_list.csv')
		new_entry = pd.DataFrame([self.log])
		log = pd.concat([log, new_entry], ignore_index=True)
		log.to_csv('cal_lists/calibrated_image_list.csv',index=False)



	def _flat_mask(self,lowlim=0.8,highlim=1.2,buffer=3):
		mask = deepcopy(self.flat) / np.nanmedian(deepcopy(self.flat))
		mask = (mask < lowlim) | (mask > highlim)
		mask = convolve(mask,np.ones((buffer,buffer)))
		mask = mask.astype(int)
		return mask

	def _saturaton_mask(self,satlimit=6e4,buffer=3):
		mask = deepcopy(self.image)
		mask = mask > satlimit
		mask = convolve(mask,np.ones((buffer,buffer)))
		mask = mask.astype(int)
		return mask


	def Make_mask(self):
		flat_mask = self._flat_mask() * 4
		sat_mask = self._saturaton_mask() * 2

		self.mask = flat_mask | sat_mask

		self._update_header_mask_bits()


	def _update_header_mask_bits(self):
		#head['STARBIT']  = (1, 'bit value for normal sources')
		self.header['SATBIT']   = (2, 'bit value for saturated sources')
		self.header['FLATBIT'] = (4, 'bit value for bad flat')
		#head['STRAPBIT'] = (8, 'bit value for bad pixels')
		#head['USERBIT']  = (16, 'bit value for USER list')
		#head['SNBIT']    = (32, 'bit value for SN list')

	
