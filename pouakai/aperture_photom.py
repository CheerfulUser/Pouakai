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
from photutils.datasets import make_100gaussians_image
from calibrimbore import sauron
from calibrimbore import get_ps1
from copy import deepcopy

import os
package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'

class ap_photom():

	def __init__(self,file=None,data=None,wcs=None,header=None,
				 fwhm=5.0,threshold=5.0,run=True,cal_model='ckmodel'):
		self.file = file
		self.data = data
		self.wcs = wcs
		self.header = header
		self.hdu = None
		self.band = None

		

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
			self.calculate_zp(fwhm,threshold)




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


	def find_sources(self,fwhm = 5, threshold=10):

		daofind = DAOStarFinder(fwhm=fwhm, threshold= threshold * self.data_std)	
		sources = daofind(self.data - self.data_median)
		print(len(sources))
		self.sources = sources


	def _calc_radii(self):
		xcoords = self.sources['xcentroid'] + 0.5
		xcoords = xcoords.astype(int)
		self.source_x = xcoords

		ycoords = self.sources['ycentroid'] + 0.5
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
				fwhm = 5 # dummy number 
			radii += [fwhm]
		self.radii = np.array(radii)


	def _get_apertures(self):
		
		xcoords = self.sources['xcentroid'] + 0.5
		xcoords = xcoords.astype(int)
		ycoords = self.sources['ycentroid'] + 0.5
		ycoords = ycoords.astype(int)
		apertures = []
		annuli = []
		pos_index = 0

		for radius in self.radii:
			position = [xcoords[pos_index],ycoords[pos_index]]
			aperture = CircularAperture(position, r=radius*1.5)
			apertures.append(aperture)
			annulus_aperture = CircularAnnulus(position, r_in=radius*3.5, r_out=radius*5)
			annuli.append(annulus_aperture)
			pos_index += 1

		self.apertures = apertures
		self.sky_ap = annuli


	def ap_photometry(self):


		masks = []
		for annulus in self.sky_ap:
			annulus_mask = annulus.to_mask(method='center')
			masks.append(annulus_mask)
			
		annulus_data = masks[0].multiply(self.data)
		
		index = masks[0].data
		annulus_data_1d = annulus_data[index > 0]
		annulus_data_1d.shape
		
		bkg_median = []
		bkg_std = []
		for mask in masks:
			annulus_data = mask.multiply(self.data)
			annulus_data_1d = annulus_data[mask.data > 0]
			mode, med, std = sigma_clipped_stats(annulus_data_1d)
			bkg_median += [med]
			bkg_std += [std]
		
		phots = []
		for i in self.apertures:
			phot = aperture_photometry(self.data, i)
			phots.append(phot)
		i = 0
		for phot in phots:
			phot['annulus_median'] = bkg_median[i]
			phot['aper_bkg'] = bkg_median[i] * self.apertures[i].area
			phot['counts'] = phot['aperture_sum'] - phot['aper_bkg']
			phot['e_counts'] = bkg_std[i] * self.apertures[i].area
			i += 1

		for phot in phots:
			for col in phot.colnames:
				phot[col].info.format = '%.8g'
			phot['sysmag'] = -2.5*np.log10(phot['counts'])
			
		table = phots[0]
		table = table.to_pandas()

		for phot in phots[1:]:
			phot = phot.to_pandas()
			table = table.append(phot, ignore_index=True)

		self.ap_photom = table


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
		ra, dec = self.wcs.all_pix2world(self.ap_photom['xcenter'],self.ap_photom['ycenter'],0)
		if self.sauron is not None:
			self.pred_mag = self.sauron.estimate_mag(ra=ra,dec=dec)

	def calc_zp(self,maglim=17,brightlim=14):
		zps = self.pred_mag - self.ap_photom['sysmag'].values
		ind = (self.pred_mag > brightlim) & (self.pred_mag < maglim)
		# cut out saturated and faint sources
		zps[~ind] = np.nan
		ind = sigma_clip(zps).mask
		zps[ind] = np.nan
		zp = np.nanmedian(zps)
		zp_std = np.nanstd(zps)
		self.zps = zps
		self.zp = zp
		self.zp_std = zp_std

	def calculate_zp(self,fwhm,threshold=10):
		self._load_image()
		self._image_stats()
		self.find_sources(fwhm=fwhm,threshold=threshold)
		self._calc_radii()
		self._get_apertures()
		self.ap_photometry()
		self._load_sauron()
		self.predict_mags()
		self.calc_zp()

	def magnitude_limit(self):
		"""Returns the magnitude limit of the filter at a given signal to noise raio"""
		sig_noise = (self.ap_photom['aper_sum_bkgsub'] / self.ap_photom['aperture_std']).values
		mag = (self.ap_photom['mag'] + self.zps).values

		ind = np.isfinite(mag) & np.isfinite(np.log10(sig_noise))
		self.snr_model =  np.polyfit(np.log10(counts_to_noise)[ind], mag[ind], 1)
		sigclip = ~sigma_clip(mag[ind] - self.fitted_line(sig_noise[ind])).mask
		fitted_model, cov = np.polyfit(np.log10(sig_noise)[ind][sigclip], mag[ind][sigclip], 1, cov=True)
		
		self.snr_model = fitted_model
		self.maglim5 = self.fitted_line(5)
		self.maglim3 = fitted_line(3)
		#chi_squared = np.sum((np.polyval(pfit, mag[ind][sigclip]) - np.log10(sig_noise[ind][sigclip])) ** 2)

		
	def mag_limit_fig(self,ax):
		sig_noise = (self.ap_photom['aper_sum_bkgsub'] / self.ap_photom['aperture_std']).values
		mag = (self.ap_photom['mag'] + self.zps).values
		ind = np.isfinite(mag) & np.isfinite(np.log10(sig_noise))
		yz = np.linspace(1,10**5,295)

		ax.plot(mag[ind],np.log10(sig_noise[ind]),'.',alpha=0.5)
		ax.plot(self.fitted_line(yz),np.log10(yz),'-')

		ax.axhline(np.log10(3),ls='-.')
		ax.axhline(np.log10(5),ls='--')
		ax.set_ylabel('log10(SNR)')
		ax.set_xlabel('Magnitude Limit')

		ax.text(15.4,3,r'$3\sigma=$ {:.2f}'.format(self.fitted_line(3)))
		ax.text(15.4,5,r'$5\sigma=$ {:.2f}'.format(self.fitted_line(5)))
	 	
		ax.set_ylim(0,10)
		ax.set_xlim(12,23)


	def fitted_line(self, sn):
		return self.snr_model[1] + self.snr_model[0] * np.log10(sn)
	
