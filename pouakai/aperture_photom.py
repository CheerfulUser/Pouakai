import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.coordinates import SkyCoord
from astropy.visualization import SqrtStretch, simple_norm, ImageNormalize
from astropy.wcs import WCS
from photutils import DAOStarFinder, aperture_photometry
from photutils.aperture import SkyCircularAperture, CircularAperture, CircularAnnulus, aperture_photometry
from photutils.aperture import ApertureStats

from scipy.optimize import minimize

from calibrimbore import sauron

from copy import deepcopy

import os
package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'

class ap_photom():

	def __init__(self,file=None,data=None,wcs=None,mask=None,header=None,ax=None,
				 threshold=5.0,run=True,cal_model='ckmodel',brightlim=14):
		self.file = file
		self.data = data
		self.wcs = wcs
		self.header = header
		self.mask = mask
		self.hdu = None
		self.band = None
		self.brightlim = brightlim


		

		self.data_median = None 
		self.data_std = None

		self.sources = None 
		self.radii = None
		self.apertures = None
		self.sky_ap = None
		self.ap_photom = None

		self.cal_sys = None
		self.cal_model = cal_model.lower()
		self.band = None
		self.sauron = None
		self.zp = None 
		self.zp_std = None  
		self.zps = None

		if run:
			self._load_image()
			self.calculate_zp(threshold)
			self.mag_limit_fig(ax)




	def _load_image(self):
		if self.file is not None:
			self.hdu = fits.open(self.file)[0]
			self.header = self.hdu.header
			self.data = self.hdu.data
			self.wcs = WCS(self.header)
			self._get_filter()
		else:
			#self._check_input()
			self._get_filter()

	def _get_filter(self):
		self.band = self.header['COLOUR'].strip(' ')



	def _image_stats(self,sigma=3):
		mode,median,std = sigma_clipped_stats(self.data, sigma=3.0)
		self.data_median = median 
		self.data_std = std


	def find_sources(self,fwhm = 3, threshold=10):

		daofind = DAOStarFinder(fwhm=fwhm, threshold= threshold * self.data_std)	
		sources = daofind(self.data - self.data_median)
		print(len(sources))
		self.sources = sources


	def _calc_radii(self):
		xcoords = self.sources['xcentroid']# + 0.5
		xcoords = xcoords.astype(int)
		self.source_x = xcoords

		ycoords = self.sources['ycentroid'] #+ 0.5
		ycoords = ycoords.astype(int)
		self.source_y = ycoords
		
		data = deepcopy(self.data) - self.data_median

		data[data < 0] = 0

		radii = []
		i = 0
		for index in range(len(xcoords)):
			data1 = (data)[ycoords[index]][xcoords[index]:xcoords[index]+100]
			normal = data1 / (data)[ycoords[index]][xcoords[index]]
			try:
				fwhm = np.where(normal < 0.5)[0][0]
			except:
				fwhm = np.nan # dummy number 
			radii += [fwhm]
		self.radii = np.array(radii) * 1.5
		self.radius = np.nanmedian(self.radii)


	def _get_apertures(self):
		
		xcoords = self.sources['xcentroid'] #+ 0.5
		ycoords = self.sources['ycentroid'] #+ 0.5
		positions = []
		for i in range(len(xcoords)):
			positions += [[xcoords[i],ycoords[i]]]
		positions = np.array(positions)
		self.aperture = CircularAperture(positions, r=self.radius)
		self.sky_ap = CircularAnnulus(positions, r_in=self.radius*3.5, r_out=self.radius*6)
		


	def ap_photometry(self):


		phot_table = aperture_photometry(self.data,self.aperture)
		aperstats = ApertureStats(self.data, self.sky_ap)
		bkg_median = aperstats.median
		bkg_std = aperstats.std
		area = self.aperture.area_overlap(self.data)
		phot_table['bkg_median'] = bkg_median
		phot_table['bkg_std'] = bkg_std
		phot_table['area'] = area
		phot_table['aper_bkg'] = area * bkg_median
		phot_table['counts'] = phot_table['aperture_sum'] - phot_table['aper_bkg']
		phot_table['e_counts'] = bkg_std * area
		phot_table['sysmag'] = -2.5*np.log10(phot_table['counts'])


		for col in phot_table.colnames:
			phot_table[col].info.format = '%.8g'  # for consistent table output

			table = phot_table.to_pandas()

		self.ap_photom = table
		near_mask = self._check_mask()
		#near_source = self._check_distance()
		ind =  (near_mask<100) #& (near_source==0)
		self.ap_photom['flag'] = ind * 0


	def _load_sauron(self):
		ra, dec = self.wcs.all_pix2world(self.ap_photom['xcenter'],self.ap_photom['ycenter'],0)
		if (dec<-35).any():
			self.cal_sys = 'skymapper'
		else:
			self.cal_sys = 'ps1'
		fname = 'cal_files/MOA-{filt}_{sys}_{model}.npy'.format(filt=self.band,
										sys=self.cal_sys,model=self.cal_model)
		self.sauron = sauron(load_state = package_directory + fname)

	def predict_mags(self):
		ind = self.ap_photom['flag'].values == 0
		ra, dec = self.wcs.all_pix2world(self.ap_photom['xcenter'].iloc[ind],self.ap_photom['ycenter'].iloc[ind],0)
		if self.sauron is not None:
			print(ra)
			mags = self.sauron.estimate_mag(ra=ra,dec=dec,close=True)
			self.pred_mag = np.ones(len(ind)) * np.nan
			self.pred_mag[ind] = mags


	def calc_zp(self,snr_lim=10):
		zps = self.pred_mag - self.ap_photom['sysmag'].values
		snr = self.ap_photom['counts'].values / self.ap_photom['e_counts'].values
		ind = (self.pred_mag > self.brightlim) & (snr > snr_lim) & (self.pred_mag < 19)
		# cut out saturated and faint sources
		zps[~ind] = np.nan
		ind = sigma_clip(zps).mask
		zps[ind] = np.nan
		zp = np.nanmedian(zps)
		zp_std = np.nanstd(zps)
		self.zps = zps
		self.zp = zp
		self.zp_std = zp_std

	def calculate_zp(self,threshold=10):
		self._load_image()
		self._image_stats()
		self.radius = 5*1.2
		for i in range(2):
			self.find_sources(fwhm=self.radius/1.5,threshold=threshold)
			self._calc_radii()
		self._get_apertures()
		self.ap_photometry()
		self._load_sauron()
		self.predict_mags()
		self.calc_zp(snr_lim=threshold)

		self.magnitude_limit(snr_lim=threshold)


	def magnitude_limit(self,snr_lim=10):
		"""Returns the magnitude limit of the filter at a given signal to noise raio"""
		sig_noise = (self.ap_photom['counts'] / self.ap_photom['e_counts']).values
		mag = (self.ap_photom['sysmag'] + self.zp).values

		ind1 = np.isfinite(mag) & np.isfinite(np.log10(sig_noise))
		ind2 = (self.pred_mag > self.brightlim) & (sig_noise > snr_lim)
		ind = ind1 & ind2
		print(len(mag[ind]))
		self.snr_model =  minimize(self._maglim_minimizer,[-1,15],args=(sig_noise[ind],mag[ind])).x
		sigclip = ~sigma_clip(mag[ind] - self.fitted_line(sig_noise[ind])).mask
		self.snr_model =  minimize(self._maglim_minimizer,[self.snr_model[0],self.snr_model[1]],
								   args=(sig_noise[ind][sigclip],mag[ind][sigclip])).x
		

		self.maglim5 = self.fitted_line(5)
		self.maglim3 = self.fitted_line(3)
		#chi_squared = np.sum((np.polyval(pfit, mag[ind][sigclip]) - np.log10(sig_noise[ind][sigclip])) ** 2)

		
	def mag_limit_fig(self,ax):
		sig_noise = (self.ap_photom['counts'] / self.ap_photom['e_counts']).values
		mag = (self.ap_photom['sysmag'] + self.zp).values
		ind = np.isfinite(mag) & np.isfinite(np.log10(sig_noise))
		sigclip = ~sigma_clip(mag[ind] - self.fitted_line(sig_noise[ind])).mask
		yz = np.linspace(1,10**5,295)

		ax.plot(mag[ind][sigclip],np.log10(sig_noise[ind][sigclip]),'.',alpha=0.5)
		ax.plot(self.fitted_line(yz),np.log10(yz),'-')

		ax.axhline(np.log10(3),ls='-.',color='k')
		ax.axhline(np.log10(5),ls='--',color='k')
		ax.set_ylabel(r'log$_{10}$(SNR)')
		ax.set_xlabel('Magnitude Limit')

		ax.text(19,2,r'$3\sigma=$ {:.2f}'.format(self.fitted_line(3)))
		ax.text(19,2.5,r'$5\sigma=$ {:.2f}'.format(self.fitted_line(5)))
	 	
		ax.set_ylim(0,3)
		ax.set_xlim(13,21)


	def _maglim_minimizer(self,var,snr,mag):
		self.snr_model = var
		mod_mag = self.fitted_line(snr)
		diff = (mag - mod_mag)**2
		return np.nansum(diff)

	def fitted_line(self, sn):
		return self.snr_model[1] + self.snr_model[0] * np.log10(sn)

	def _check_mask(self):
		if self.mask is not None:
			flags = np.where(self.mask > 0)
			bx = flags[1]
			by = flags[0]

			sx = self.ap_photom['xcenter'].values
			sy = self.ap_photom['ycenter'].values
			r = self.radii 

			x = (sx[:,np.newaxis] - bx[np.newaxis,:])
			y = (sy[:,np.newaxis] - by[np.newaxis,:])

			d2 = x**2 + y**2

			near = d2 < r[:,np.newaxis]
			near = np.nansum(near,axis=1)
		else:
			near = self.radii * 0
		return near.astype(int)

	def _check_distance(self,limit=20):

		sx = self.ap_photom['xcenter'].values
		sy = self.ap_photom['ycenter'].values

		x = sx[:,np.newaxis] - sx[np.newaxis,:]
		y = sy[:,np.newaxis] - sy[np.newaxis,:]

		d2 = x**2 + y**2
		d2[d2==0] = np.nan

		mins = np.nanmin(d2,axis=1)
		mins = mins < limit**2

		return mins

