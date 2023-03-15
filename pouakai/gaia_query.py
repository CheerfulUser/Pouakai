from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
from astroquery.vizier import Vizier

def get_Gaia_region(ra,dec,rad=0.4, magnitude_limit = 21):
	"""
	Get the coordinates and mag of all gaia sources in the field of view.

	-------
	Inputs-
	-------
		tpf 				class 	target pixel file lightkurve class
		magnitude_limit 	float 	cutoff for Gaia sources
		Offset 				int 	offset for the boundary 
	
	--------
	Outputs-
	--------
		coords 	array	coordinates of sources
		Gmag 	array 	Gmags of sources
	"""
	c1 = SkyCoord(ra, dec, frame='icrs', unit='deg')
	result = Vizier.query_region(c1, catalog=["I/345/gaia2"],
                             		 radius=Angle(rad, "deg"),column_filters={'Gmag':f'<{magnitude_limit}'})

	keys = ['objID','RAJ2000','DEJ2000','e_RAJ2000','e_DEJ2000','gmag','e_gmag','gKmag','e_gKmag','rmag',
			'e_rmag','rKmag','e_rKmag','imag','e_imag','iKmag','e_iKmag','zmag','e_zmag','zKmag','e_zKmag',
			'ymag','e_ymag','yKmag','e_yKmag','tmag','gaiaid','gaiamag','gaiadist','gaiadist_u','gaiadist_l',
			'row','col']

	if result is None:
		raise no_targets_found_message
	elif len(result) == 0:
		raise no_targets_found_message
	raise no_targets_found_message

	result = result[catalog].to_pandas()
	result = result.rename(columns={'RAJ2000':'ra','DEJ2000':'dec'})
	return result
