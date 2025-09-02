from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
import pandas as pd
from astroquery.vizier import Vizier

def get_gaia_region(ra,dec,size=0.4, magnitude_limit = 21):
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
	c1 = SkyCoord(ra, dec, unit='deg')
	Vizier.ROW_LIMIT = -1

	result = Vizier.query_region(c1, catalog=["I/345/gaia2"],
                             		 radius=Angle(size, "arcsec"),column_filters={'Gmag':f'<{magnitude_limit}'})

	keys = ['objID','RAJ2000','DEJ2000','e_RAJ2000','e_DEJ2000','gmag','e_gmag','gKmag','e_gKmag','rmag',
			'e_rmag','rKmag','e_rKmag','imag','e_imag','iKmag','e_iKmag','zmag','e_zmag','zKmag','e_zKmag',
			'ymag','e_ymag','yKmag','e_yKmag','tmag','gaiaid','gaiamag','gaiadist','gaiadist_u','gaiadist_l',
			'row','col']


	no_targets_found_message = ValueError('Either no sources were found in the query region '
                                          'or Vizier is unavailable')
	if result is None:
		raise no_targets_found_message
	elif len(result) == 0:
		raise no_targets_found_message
	

	result = result['I/345/gaia2'].to_pandas()
	#result = result.rename(columns={'RA_ICRS':'ra','DE_ICRS':'dec'})
	#account for proper motion
	return result

