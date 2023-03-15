import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.coordinates import SkyCoord
from astropy.visualization import SqrtStretch, simple_norm, ImageNormalize
from astropy.wcs import WCS
from photutils import DAOStarFinder, aperture_photometry
from photutils.aperture import SkyCircularAperture, CircularAperture, CircularAnnulus, aperture_photometry, ApertureStats

from scipy.ndimage.filters import convolve

from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata

from mpl_toolkits.mplot3d import Axes3D


from scipy.optimize import minimize

from calibrimbore import sauron, get_skymapper_region, get_ps1_region

from copy import deepcopy
from .gaia_query import get_gaia_region

import os
package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'

class cal_photom():

	def __init__(self,file=None,data=None,wcs=None,mask=None,header=None,ax=None,
				 threshold=5.0,run=True,cal_model='ckmodel',brightlim=10,rescale=True,
				 plot=True,floor=None,radius_override=None,use_catalogue=True,
				 band_override=None, ):
		self.file = file
		self.data = data
		self.wcs = wcs
		self.header = header
		self.mask = mask
		self.hdu = None
		self.band = None
		self.brightlim = brightlim
		self.image_floor = floor
		self.radius_override = radius_override
		self.use_catalogue = use_catalogue
		self._band_override = band_override


		

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
		self.phot_table = None

		if run:
			self._load_image()
			self._basic_source_mask()
			self.calculate_zp(threshold)
			if rescale:
				self.ZP_correction()
				self.Recast_image_scale()
				self.calculate_zp(threshold)
			else:
				self.zp_surface = np.ones_like(self.data)
			self.ap_photom['mag'] = self.ap_photom['sysmag'] + self.zp
			self.ap_photom['e_mag'] = 2.5/np.log(10)*(self.ap_photom['e_counts']/self.ap_photom['counts'])
			if plot:
				self.mag_limit_fig(ax)

			self.ap_photom['gaiaID'] = self.cat['objID'].values




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
		if self._band_override is None:
			self.band = self.header['FILTER'].strip(' ')
			if self.band in {'B','V','R','I'}:
				self.band = 'bessell_' + self.band
			if self.band in {'g','r','i','z'}:
				self.band = 'sloan_' + self.band
		else:
			self.band = self._band_override



	def _image_stats(self,sigma=3):
		mode,median,std = sigma_clipped_stats(self.data, sigma=sigma)
		self.data_median = median 
		self.data_std = std


	def find_sources(self,fwhm = 5, threshold=100):

		daofind = DAOStarFinder(fwhm=fwhm, threshold= threshold * self.data_std)	
		sources = daofind(self.data - self.data_median)
		self.sources = sources.to_pandas()
		self._mask_killer()

	def catalogue_sources(self):
		ra,dec = self.wcs.all_pix2world(self.data.shape[1]//2,self.data.shape[0]//2,0)
		#if dec > -25:
		#	cat = get_ps1_region([ra],[dec],size=.4*60**2)
		#else:
		#	cat = get_skymapper_region([ra],[dec],size=.4*60**2)
		cat = get_gaia_region([ra],[dec],size=.2*60**2)
		tab = deepcopy(cat)#.iloc[np.isfinite(cat.r.values)]
		x, y = self.wcs.all_world2pix(tab.ra.values,tab.dec.values,0)
		tab['x'] = x
		tab['y'] = y

		ind = (x > 15) & (x < self.data.shape[1]-15) & (y > 15) & (y < self.data.shape[0]-15)
		tab = tab.iloc[ind]
		self.cat = tab
		sources = tab[['x','y']]
		sources.rename(columns={'x': 'xcentroid', 'y': 'ycentroid'}, inplace=True)
		self.sources = sources
		self._mask_killer()


	def _calc_radii(self):
		xcoords = self.sources['xcentroid'].values + 0.5
		self.source_x = xcoords
		xcoords = xcoords.astype(int)
		
		ycoords = self.sources['ycentroid'].values + 0.5
		self.source_y = ycoords
		ycoords = ycoords.astype(int)
		
		data = deepcopy(self.data) - self.data_median

		data[data < 0] = 0

		radii = []

		for i in range(len(xcoords)):
			#print(xcoords[i],ycoords[i])
			data1 = data[ycoords[i],xcoords[i]:xcoords[i]+30]
			normal = data1 / data[ycoords[i], xcoords[i]]
			try:
				fwhm = np.where(normal < 0.5)[0][0]
			except:
				fwhm = np.nan # dummy number 
			radii += [fwhm]
		self.radii = np.array(radii) * 1.4
		self.radius = np.nanmedian(self.radii)


	def _get_apertures(self):
		
		xcoords = self.sources['xcentroid'].values
		ycoords = self.sources['ycentroid'].values
		positions = []
		for i in range(len(xcoords)):
			positions += [[xcoords[i],ycoords[i]]]
		positions = np.array(positions)
		if self.radius_override is not None:
			self.aperture = CircularAperture(positions, r=self.radius_override)
		else:
			self.aperture = CircularAperture(positions, r=self.radius)
		self.sky_ap = CircularAnnulus(positions, r_in=self.radius*2, r_out=self.radius*10)
		


	def ap_photometry(self):
		phot_table = aperture_photometry(self.data,self.aperture)
		masked = deepcopy(self.data)
		masked[self.source_mask==0] = np.nan
		aperstats = ApertureStats(masked, self.sky_ap)
		unmasked = ApertureStats(self.data, self.sky_ap)
		bkg_median = aperstats.median
		bkg_std = aperstats.std
		ind = np.isnan(bkg_median)
		bkg_median[ind] = unmasked.median[ind]
		bkg_std[ind] = unmasked.std[ind]
		area = self.aperture.area_overlap(self.data)
		phot_table['bkg_median'] = bkg_median
		phot_table['bkg_std'] = bkg_std 
		phot_table['area'] = area
		phot_table['aper_bkg'] = area * bkg_median
		phot_table['counts'] = phot_table['aperture_sum'] - phot_table['aper_bkg']
		phot_table['e_counts'] = phot_table['bkg_std'] * area
		phot_table['snr'] = phot_table['counts'] / phot_table['e_counts']
		phot_table['sysmag'] = -2.5*np.log10(phot_table['counts'])


		for col in phot_table.colnames:
			phot_table[col].info.format = '%.8g'  # for consistent table output

			table = phot_table.to_pandas()

		self.ap_photom = table
		#near_mask = self._check_mask()
		#near_source = self._check_distance()
		#ind =  (near_mask<100) #& (near_source==0)
		self.ap_photom['flag'] = 0 #ind


	def _load_sauron(self):
		ra, dec = self.wcs.all_pix2world(self.ap_photom['xcenter'],self.ap_photom['ycenter'],0)
		if (dec < -25).any():
			self.cal_sys = 'skymapper'
		else:
			self.cal_sys = 'ps1'
		fname = 'cal_files/{filt}_{sys}_{model}.npy'.format(filt=self.band,
										sys=self.cal_sys,model=self.cal_model)
		self.sauron = sauron(load_state = package_directory + fname)

	def predict_mags(self):
		ind = self.ap_photom['flag'].values == 0
		ra, dec = self.wcs.all_pix2world(self.ap_photom['xcenter'].iloc[ind],self.ap_photom['ycenter'].iloc[ind],0)
		if self.sauron is not None:
			mags = self.sauron.estimate_mag(ra=ra,dec=dec,close=True)
			self.pred_mag = np.ones(len(ind)) * np.nan
			self.pred_mag[ind] = mags


	def calc_zp(self,snr_lim=3):
		zps = self.pred_mag - self.ap_photom['sysmag'].values
		snr = self.ap_photom['snr'].values
		ind = (self.pred_mag > self.brightlim) & (snr > snr_lim) & (self.pred_mag < 17)
		# cut out saturated and faint sources
		zps[~ind] = np.nan
		ind = sigma_clip(zps,maxiters=10).mask
		zps[ind] = np.nan
		zp = np.nanmedian(zps)
		zp_std = np.nanstd(zps)
		self.zps = zps
		self.zp = zp
		self.zp_std = zp_std
		self.good = ind

	def calculate_zp(self,threshold=10):
		self._load_image()
		self._image_stats()
		
		if self.use_catalogue:
			self.catalogue_sources()
			self._calc_radii()
		else:
			self.radius = 3*1.2
			for i in range(2):
				self.find_sources(fwhm=self.radius/1.2,threshold=threshold)
				self._calc_radii()

		self._get_apertures()
		self.ap_photometry()
		self._load_sauron()
		self.predict_mags()
		self.calc_zp(snr_lim=threshold)

		self.magnitude_limit(snr_lim=threshold)


	def magnitude_limit(self,snr_lim=10):
		"""Returns the magnitude limit of the filter at a given signal to noise raio"""
		sig_noise = self.ap_photom['snr'].values
		mag = (self.ap_photom['sysmag'] + self.zps).values

		ind1 = np.isfinite(mag) & np.isfinite(np.log10(sig_noise))
		ind2 = (self.pred_mag > self.brightlim) & (sig_noise > snr_lim)
		ind = ind1 & ind2
		self.snr_model =  minimize(self._maglim_minimizer,[-1,15],args=(sig_noise[ind],mag[ind])).x
		sigclip = ~sigma_clip(mag[ind] - self.fitted_line(sig_noise[ind])).mask
		self.snr_model =  minimize(self._maglim_minimizer,[self.snr_model[0],self.snr_model[1]],
								   args=(sig_noise[ind][sigclip],mag[ind][sigclip])).x
		

		self.maglim5 = self.fitted_line(5)
		self.maglim3 = self.fitted_line(3)
		#chi_squared = np.sum((np.polyval(pfit, mag[ind][sigclip]) - np.log10(sig_noise[ind][sigclip])) ** 2)

		
	def mag_limit_fig(self,ax):
		sig_noise = self.ap_photom['snr'].values
		mag = (self.ap_photom['sysmag'] + self.zps).values
		ind = np.isfinite(mag) & np.isfinite(np.log10(sig_noise))
		sigclip = ~sigma_clip(mag[ind] - self.fitted_line(sig_noise[ind])).mask
		yz = np.linspace(1,10**5,295)

		#ax.plot(mag[ind],np.log10(sig_noise[ind]),'.',alpha=0.5)
		if ax is not None:
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

	def _mask_killer(self,buffer=5):
		if self.mask is not None:
			mask = convolve(self.mask,np.ones((buffer,buffer)))
			x,y = np.where(mask > 0)
			sx = self.sources['xcentroid'].values
			sy = self.sources['ycentroid'].values
			d = np.sqrt((sx[:,np.newaxis] - x[np.newaxis,])**2 + (sy[:,np.newaxis] - y[np.newaxis,])**2)
			mind = np.nanmin(d,axis=1)
			ind = mind > 1
			print(f'Killed {len(sx) - np.sum(ind*1)} sources')
			self.sources = self.sources.iloc[ind]
			
			self.source_x = self.sources['xcentroid'].values
			self.source_y = self.sources['ycentroid'].values

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

	def _basic_source_mask(self,buffer=0.04):
		"""
		Basic method to mask pixels that are significantly elevated from the background.
		Most useful for not including stars in background determination.
		1 is background.
		"""
		med = np.nanmedian(self.data)
		limit = med + (buffer * med)

		mask = self.data < limit
		self.source_mask = mask

	def Fit_surface(self,mask=None,smoother=100):
		ind = np.isfinite(self.zps)
		if mask is not None:
			ind = ind & mask
		x_data = (self.ap_photom['xcenter'].values[ind] + 0.5).astype(int)
		y_data = (self.ap_photom['ycenter'].values[ind] + 0.5).astype(int)
		z_data = self.zps[ind]
		
		zpimage = np.zeros_like(self.data)
		zpimage[y_data,x_data] = z_data
		zpimage[zpimage==0] = np.nan
		#zpimage[cut.mask] = np.nan

		x = np.arange(0, zpimage.shape[1])
		y = np.arange(0, zpimage.shape[0])
		arr = np.ma.masked_invalid(zpimage)
		xx, yy = np.meshgrid(x, y)
		#get only the valid values
		x1 = xx[~arr.mask]
		y1 = yy[~arr.mask]
		newarr = arr[~arr.mask]

		estimate = griddata((x1, y1), newarr.ravel(),
									(xx, yy),method='linear')
		bitmask = np.zeros_like(zpimage,dtype=int)
		bitmask[np.isnan(estimate)] = 128 | 4
		nearest = griddata((x1, y1), newarr.ravel(),
									(xx, yy),method='nearest')

		estimate[np.isnan(estimate)] = nearest[np.isnan(estimate)]

		estimate = gaussian_filter(estimate,smoother)

		return estimate, bitmask

	def Recast_image_scale(self,newzp=25):
		"""
		Recast the image scale to a new zeropoint.
		"""
		if self.image_floor is not None:
			new_image = ((self.data - np.nanmedian(self.image_floor)) * 10**((self.zp_surface - newzp) / -2.5)) + np.nanmedian(self.image_floor)
		else:
			new_image = ((self.data - np.nanmedian(self.data)) * 10**((self.zp_surface - newzp) / -2.5)) + np.nanmedian(self.data)
		self.data = new_image


	def ZP_correction(self,sigma=2):
		"""
		Correct the zeropoint for the residual background varibility.
		"""
		#tmp,_ = self.Fit_surface(mask=None,smoother=200)
		#x_data = (self.ap_photom['xcenter'].values + 0.5).astype(int)
		#y_data = (self.ap_photom['ycenter'].values + 0.5).astype(int)
		#diff = (self.zps - tmp[y_data.astype(int),x_data.astype(int)])
		#cut = ~sigma_clip(diff,sigma=sigma).mask
		#estimate,bitmask = self.Fit_surface(mask=cut,smoother=30)
		self.zp_surface = np.nanmedian(self.zps)
		


	def plot_zp_correction(self,ax,cut):
		"""
		Deosn't seem likek this will work without some restructuring.
		Plots the zeropoint correction for the image.
		"""
		x_data = (self.ap_photom['xcenter'].values + 0.5).astype(int)
		y_data = (self.ap_photom['ycenter'].values + 0.5).astype(int)
		x = np.arange(0, self.data.shape[1])
		y = np.arange(0, self.data.shape[0])
		xx, yy = np.meshgrid(x, y)
		ax = Axes3D(ax)
		# plot surface
		#ax.plot_surface(xx, yy, tmp,alpha=0.5,color='C1')
		surf = ax.plot_surface(xx, yy, self.zp_surface,alpha=0.3,color='C1',label='Correction surface')
		surf._facecolors2d = surf._facecolor3d
		surf._edgecolors2d = surf._edgecolor3d
		# plot input data
		ax.scatter(x_data[~cut], y_data[~cut], self.zps[~cut], color='C3',label='Outliers')
		ax.scatter(x_data[cut], y_data[cut], self.zps[cut], color='C0',label='Calibrations')
		# set plot descriptions
		ax.set_xlabel('x',fontsize=15)
		ax.set_ylabel('y ',fontsize=15)
		ax.set_zlabel('Zeropoint',fontsize=15)
		ax.legend()
		ax.view_init(elev=13.,azim=-167)
		plt.tight_layout()