from astropy.io import fits
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from astroquery.astrometry_net import AstrometryNet

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


		if (reduction.lower() == 'full'):
			self.reduce_image()

			self.save_intermediate()


	def _read_fits(self):
		hdu = fits.open(file)[0]
		self.header = hdu.header
		self.raw_image = hdu.data
		self.jd = hdu.header['JDSTART']
		self.filter = hdu.header['COLOUR']
		self.chip = hdu.header['CHIP']
		self.exp_time = hdu.header['EXPTIME']

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
		hdu = fits.open(file)[0]
		data = hdu.data

		if cal_type.lower() == 'flat':
			self.flat_file = file 
			self.flat = data
		elif cal_type.lower() == 'dark':
			self.dark_file = file 
			self.dark = data


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
	    

	def _get_wcs(self):
	    attempt = 0
	    solved = False
	    while (attempt < 10) | (not solved):
	        try:
	            wcs_head = ast.solve_from_image(self.tmp_file)
	            del wcs_head['COMMENT']
	            del wcs_head['HISTORY']
	            solved = True
	        except:
	            attempt += 1
	    if (attempt > 10) | (not solved):
	        raise ValueError('Could not solve WCS in {} attempts'.format(attempt))
	    new_head = deepcopy(header)
	    for key in wcs_head:
	        new_head[key] = (wcs_head[key],wcs_head.comments[key])
	    self.header = new_head
	    

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


    def reduce_image(self):
    	self._check_vars()

    	image = (self.raw_image - self.dark) / self.flat

    	self.image = image


	def save_intermediate(self):
		self._update_header_sky()
		self._update_header_dark()
		self._update_header_flat()
		name = self.savepath + 'red/' + self.base_name + '_red.fits.gz'
		self.red_name = name


		if verbose:
			print('Saving intermediated calibrated file')
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

    def save_intermediate_wcs(self):
    	name = self.savepath + 'wcs/' + self.base_name + '_wcs.fits.gz'
		self.wcs_name = name

		if verbose:
			print('Saving intermediated wcs file')
		fits.writeto(name,self.image,header=self.header,overwrite=True)




	def _get_system(self):
		low = min(self.wcs.calc_footprint()[:,1])
		if low < -30:
		    system = 'skymapper'
		else:
		    system = 'ps1'
	    return system

	def calculate_zp(self,):
		system = self._get_system()
		state_name = './cal_files/MOA-{}_{}.npy'.format(self.filter,system)
		moa = cal.sauron(load_state=state_name)











