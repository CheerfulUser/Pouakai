from astropy.io import fits
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from astroquery.astrometry_net import AstrometryNet
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.stats import iqr
from aperture_photom import ap_photom

class pouakai():

	def __init__(self,file,reduction='full',time_tolerence=1,dark_tolerence=1,savepath=''):
		self.file = file 
		self.savepath = savepath
		self._set_base_name()
		self.time_tolerence = time_tolerence
		self.dark_tolerence = dark_tolerence
		self._read_fits()

		self._get_master('dark')
		self._get_master('flat')

		self._setup_fig()

		self.reduce(reduction)

	def reduce(self,reduction,local_astrom=True):
		self._check_reduction(reduction)
		self.reduce_image()
		if local_astrom:
			self.wcs_astrometrynet_local()
		else:
			self.wcs_astrometrynet()

		self.calculate_zp()
		self.savefig()
		if (reduction == 'full'):
			self.save_intermediate()
			

	def _start_record(self):
		document = {'name':None,'band':None,'chip':None,
					'exptime':None,'jd':None,'date':None,
					'field':None,'filename':None,'flat':None,
					'dark':None,'tdiff_flat':None,
					'tdiff_dark':None,'zp':None,'zperr':None,
					'maglim5':None,'maglim3':None}


	def _read_fits(self):
		hdu = fits.open(file)[0]
		self.header = hdu.header
		self.raw_image = hdu.data
		self.jd = hdu.header['JDSTART']
		self.filter = hdu.header['COLOUR']
		self.chip = hdu.header['CHIP']
		self.exp_time = hdu.header['EXPTIME']
		self._field_coords()

	def _field_coords(self):
		ra = self.header['RA'].strip(' ')
		dec = self.header['DEC'].strip(' ')
		c = SkyCoord(ra,dec, unit=(u.hourangle, u.deg))
		self.field_coord = c



	def _set_base_name(self):
		self.base_name = self.file.split('/')[-1].split('.f')[0]

	def _get_master(self,cal_type):
		
		if cal_type.lower() == 'flat':
			masters = pd.read_csv('filename')

		elif cal_type.lower() == 'dark':
			masters = pd.read_csv('filename')
			exptimes = masters['exptime'].values
			ind = np.where(abs(self.exp_time - exptimes)<self.dark_tolerence)[0]
			if len(ind) == 0:
				m = 'No master darks with exptime {}'.format(self.exp_time)
				raise ValueError(m)
			masters = masters.iloc[ind]

		else:
			raise ValueError('Only flat and dark are valid options!!')
		file = self._find_master(masters)
		hdu = fits.open(file)
		data = hdu[0].data
		err =  hdu[1].data


		if cal_type.lower() == 'flat':
			self.flat_file = file 
			self.flat = data
			self.flat_err = err
		elif cal_type.lower() == 'dark':
			self.dark_file = file 
			self.dark = data
			self.dark_err = err


	def _find_master(self,masters):
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

		return file




	def _update_header_sky(self):
		"""
		Update the fits header to the format needed for photpipe
		"""
		mean, med, std = sigma_clipped_stats(self.image, sigma=3.0)
		self.header['SKYADU'] = (med, 'median sky')
		self.header['SKYSIG'] = (std, 'median sky noise')
			
	def _update_header_dark(self):
		self.header['DARKCORR'] = (True,self.dark_file)
			
	def _update_header_flat(self):
		self.header['FLATCORR'] = (True,self.flat_file)
		
	def _update_header_zeropoint(self):
		self.header['ZP'] = (self.zp,'Calibrimbore ' + self.system)
		
		

	def _check_vars(self):
		if self.flat is None:
			m = 'No flat specified'
			raise ValueError(m)
		if self.dark is None:
			m = 'No dark specified'
			raise ValueError(m)
		if self.raw_image is None:
			m = 'No image specified'
			raise ValueError(m)

	def _check_dirs(self):
		dirlist = ['wcs','wcs_tmp','red','cal','fig']
		for d in dirlist:
			if ~os.path.isdir(self.savepath + d):
				os.mkdir(self.savepath + d)

	def reduce_image(self):
		self._check_vars()
		
		image = (self.raw_image - self.dark) / self.flat

		self.image = image

		self._add_image(self.raw_image,'A')
		self._add_image(self.flat,'B')
		self._add_image(image,'C')
		self._add_image(image,'D')


	def save_intermediate(self):
		self._update_header_sky()
		self._update_header_dark()
		self._update_header_flat()
		name = self.savepath + 'red/' + self.base_name + '_red.fits'
		self.red_name = name


		if self.verbose:
			print('Saving intermediated calibrated file')
		fits.writeto(name,self.image,header=self.header,overwrite=True)


	def save_image(self):
		self._update_header_sky()
		self._update_header_dark()
		self._update_header_flat()
		name = self.savepath + 'cal/' + self.base_name + '_cal.fits'
		self.cal_name = name

		if self.verbose:
			print('Saving final calibrated image')
		fits.writeto(name,self.image,header=self.header,overwrite=True)

	def wcs_astrometrynet(self,timeout=120):
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
		self.header = header
		self.wcs = WCS(self.header)


	def wcs_astrometrynet_local(self):
		astrom_call = "solve-field -O -o {savename} -p --ra {ra} --dec {dec} --radius 2 {file}"

		name = self.savepath + 'wcs_tmp/' + self.base_name + '_wcs'
		solver = astrom_call.format(savename = name, ra = self.field_coord.ra.deg,
									dec = self.field_coord.dec.deg, file = self.red_name)
		os.system(solver)

		wcs_header = fits.open(name + '.new')[0].header
		del wcs_head['COMMENT']
		del wcs_head['HISTORY']
		self.header = wcs_header
		self.wcs = WCS(self.header)
		
		clear = 'rm ' + self.savepath + 'wcs_tmp/' + self.base_name + '*'
		os.system(clear)
		
	


	def save_intermediate_wcs(self):
		name = self.savepath + 'wcs/' + self.base_name + '_wcs.fits'
		self.wcs_name = name

		if verbose:
			print('Saving intermediated wcs file')
		fits.writeto(name,self.image,header=self.header,overwrite=True)

	def calculate_zp(self,threshold=100,model='ckmodel'):
		self.cal = ap_photom(data=self.data,wcs=self.wcs, header=self.header,
									threshold=threshold,cal_model=model)
		self.header['ZP'] = (self.cal.zp, 'Calibrimbore zeropoint')
		self.header['ZPERR'] = (self.cal.zp_std, 'Calibrimbore zeropoint error')
		self.header['MAGLIM5'] = (self.cal.maglim5, '5 sigma mag lim')
		self.header['MAGLIM3'] = (self.cal.maglim3, '3 sigma mag lim')
		self.fig_axis['D'].plot(self.cal.source_x,source_y,'r.')
		self._zp_hist()
		self._zp_color()

	def _zp_hist(self):
		zps = self.cal.zps
		b = int(np.nanmax(zps) - np.nanmin(zps) /(2*iqr(zps)*len(zps)**(-1/3)))
		self.fig_axis['D'].hist(zps,bins=b,alpha=0.5)		
		med = self.cal.zp
		high = self.cal.zp+self.cal.zp_std
		low = self.cal.zp-self.cal.zp_std
		self.fig_axis['D'].axvline(med,color='k',ls='--')
		self.fig_axis['D'].axvline(low,color='k',ls=':')
		self.fig_axis['D'].axvline(high,color='k',ls=':')
		s = ('$zp='+str((np.round(med,2)))+'^{+' + 
			str(np.round(high-med,2))+'}_{'+
			str(np.round(low-med,2))+'}$')
		self.fig_axis['D'].annotate(s,(.75,.8),fontsize=10,xycoords='axes fraction')
		self.fig_axis['D'].set_xlabel('zeropoint',fontsize=15)
		self.fig_axis['D'].set_ylabel('Occurrence',fontsize=15)

	def _zp_color(self):
		zps = self.cal.zps
		gr = self.cal.cat_mags['g'] - self.cal.cat_mags['r']

		self.fig_axis['F'].plot(gr,zps,'.')
		self.fig_axis['F'].set_ylabel('zeropoint',fontsize=15)
		self.fig_axis['F'].set_xlabel('$g-r$',fontsize=15)	

	def _maglim_fig(self):
		r = self.cal.cat_mags['r']


	def save_fig(self):
		name = self.savepath + 'wcs_tmp/' + self.base_name + '_diag.pdf'
		self.fig.savefig(name)

	def _setup_fig(self):
		self.fig = plt.figure(figsize=(8.27,11.69),constrained_layout=True)
		self.fig_axis = fig.subplot_mosaic(
											"""
											AB
											AB
											CD
											CD
											EF
											"""
										  )
		self.fig_axis['A'].set_title('Raw image',fontsize=15)
		self.fig_axis['B'].set_title('Flat image',fontsize=15)
		self.fig_axis['C'].set_title('Reduced image',fontsize=15)
		self.fig_axis['D'].set_title('Calibration sources',fontsize=15)
		self.fig_axis['E'].set_title('Zeropoint',fontsize=15)
		self.fig_axis['F'].set_title('Color dependance',fontsize=15)

	def _add_image(self,image,ax_ind):
		vmin = np.percentile(image,16)
		vmax = np.percentile(image,84)
		self.fig_axis[ax_ind].imshow(image,origin='lower',
									 vmin=vmin,vmax=vmax)







