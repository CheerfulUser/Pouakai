from astropy.io import fits
from datetime import datetime
from astroquery.imcce import Skybot
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import sys
import os
from astropy.coordinates import Angle


def extracts_header_info(file_path):
    """
    Extracts relevant information from the header and data of a FITS file.

    Parameters:
    file_path (str): The path to the FITS file.

    Returns:
    tuple: A tuple containing the extracted information:
           (ra, dec, mag, snr, date, target, band, observatory, mag_limit)
    """
    # Open the FITS file
    hdul = fits.open(file_path)

    # Access the header and data
    header = hdul[1].header
    data = hdul[1].data

    # Extract relevant information
    ra = header['RA']
    dec = header['DEC']
    mag_list = data['mag']
    date = header['DATE']
    band = header['FILTER'].strip()
    mag_limit = header['MAGLIM5']

    # Close the FITS file
    hdul.close()

    info = ra, dec, date.strip(), band, mag_limit, mag_list
    return info, data



def get_solar_system_objects(ra, dec, observation_time, mag_limit):
    """
    Retrieves solar system objects by querying the Minor Planet Center 
    within a specified field of view and filters them based on magnitude.
    Additionally prints them as a table.

    Parameters:
    ra (float): Right Ascension of the field center in hours.
    dec (float): Declination of the field center in degrees.
    observation_time (str): Observation time in ISO format.
    mag_limit (float): Magnitude limit for filtering objects.

    Returns:
    Table: A table of solar system objects within the field of view that meet the magnitude limit criteria.
    """
    # Convert RA and Dec to degrees
    field_center = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    
    # Create Time object for the observation time
    epoch = Time(observation_time, format='isot', location='474')
    
    # Perform the cone search with a radius of 30 arcminutes
    results = Skybot.cone_search(field_center, 120 * u.arcmin, epoch)
    
    # Convert the 'V' column to float and filter results
    results['V'] = results['V'].astype(float)
    filtered_results = results[results['V'] <= float(mag_limit) * u.mag]
    
    return filtered_results


def compare_objects(data_pouakai, data_MPC, tolerance):
    """
    Compares the source list from Pouakai with the list of objects from the MPC to find coordinate matches.

    Parameters:
    data_pouakai (list): List of objects from Pouakai, each containing RA and DEC values and additional info.
    data_MPC (list): List of objects from the MPC, each containing RA and DEC values and additional info.
    tolerance (float): The maximum allowable separation (in degrees) for objects to be considered a match.

    Returns:
    list: A list of tuples, each containing a matching pair of objects from data_pouakai and data_MPC, 
    along with the index of the matched object from data_pouakai. Returns an empty list and prints a message if no matches.

    """
    matches = []
    i = -1

    for obj1 in data_pouakai:          #Iterate and compare
        ra1, dec1 = obj1[3], obj1[4]
        i += 1
        coord1 = SkyCoord(ra=ra1, dec=dec1, unit='deg')
        for obj2 in data_MPC:
            ra2, dec2 = obj2[2], obj2[3]
            coord2 = SkyCoord(ra=ra2, dec=dec2, unit='deg')
            if coord1.separation(coord2).deg <= tolerance:
                matches.append((obj1, obj2, i))
    
    if len(matches) == 0:
        print('No matches found!')
    
    return matches





def get_object_designation(matches):
    """
    Extracts relevant information from the match list and returns the formatted object designation for a minor planet.

    Parameters:
    matches (list): A list of tuples containing matched objects from Pouakai and MPC.

    Returns:
    str: The packed provisional designation for the minor planet, including the cycle number.

    """
    # From the MPC side of the matches, extract the periodic comet number and provisional designation
    MP_number = str(matches[0][1][0])  # Convert numpy.int32 to string
    designation = matches[0][1][1]

    def pack_MP_number(MP_number):
        """
        Formats the packed minor planet number into A5.

        Parameters:
        MP_number (str): The minor planet number to be formatted.

        Returns:
        str: The packed provisional designation for the minor planet.

        NOTE: There must be a better way to do this!
        """
        # Format the minor planet number according to the rules
        if MP_number.isdigit():
            number = int(MP_number)
            if number < 100000:
                MP_number = f"{number:05d}"
            elif number <= 619999:
                mod = number % 10000
                div = number // 10000
                if 10 <= div <= 35:
                    MP_number = f"{chr(ord('A') + div - 10)}{mod:04d}"
                elif 36 <= div <= 61:
                    MP_number = f"{chr(ord('a') + div - 36)}{mod:04d}"
            else:
                number -= 620000
                base62 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
                packed_number = ""
                for _ in range(4):
                    packed_number = base62[number % 62] + packed_number
                    number //= 62
                MP_number = f"~{packed_number}"

        return MP_number


    def pack_provisional_designation(designation):
        """
        Packs a provisional designation into the Minor Planet Center's packed format.

        Parameters:
        designation (str): The provisional designation to be packed (e.g., "2007 TA418").

        Returns:
        str: The packed provisional designation (e.g., "K07Tf8A").

        PROBLEM:
        Doesn't account for Survey designations of the form 
        2040 P-L, 3138 T-1, 1010 T-2 and 4101 T-3 are packed differently. 
        Columns 1-3 contain the code indicating the survey and columns 4-7 
        contain the number within the survey.
        """
        # Define the mapping for the first two digits of the year
        year_mapping = {'18': 'I', '19': 'J', '20': 'K'}

        # Split the designation into its components
        year, half_month_order = designation.split()
        
        # Pack the year
        packed_year = year_mapping[year[:2]] + year[2:]

        def pack_number(number):       #Note: there must be a better way to do this!
                if number < 100:
                    return f"{number:02d}"
                elif number <= 619:
                    mod = number % 10
                    div = number // 10
                    if 10 <= div <= 35:
                        return f"{chr(ord('A') + div - 10)}{mod:01d}"
                    elif 36 <= div <= 61:
                        return f"{chr(ord('a') + div - 36)}{mod:01d}"
                else:
                    print("Number too high to assign letter")
                    return str(number)

        letters = ''.join(filter(str.isalpha, half_month_order))
        digits = ''.join(filter(str.isdigit, half_month_order))
        if digits:
            number = int(digits)
        else:
            number = 0
        packed_number = pack_number(number)
        if letters:
            formatted_string = f"{letters[0]}{packed_number}{letters[1:]}"
        else:
            formatted_string = f"{packed_number}"

        # Combine the components into the packed designation
        packed_designation = f"{packed_year}{formatted_string}"

        return packed_designation

    # Combine the information into the correct order for the packed provisional designation
    packed_MP_number = pack_MP_number(MP_number)
    packed_provisional_designation = pack_provisional_designation(designation)

    return packed_MP_number + packed_provisional_designation


def convert_utc_format(date_str):
    # Parse the input date string
    dt = datetime.strptime(date_str, '%y-%m-%dT%H:%M:%S')
    
    # Calculate the decimal day
    decimal_day = dt.day + (dt.hour * 3600 + dt.minute * 60 + dt.second) / 86400
    
    # Format the date to the desired output with precision of 0.00001 days
    formatted_date = dt.strftime('%Y %m %d.') + f'{decimal_day%1:.5f}'[2:]
    return formatted_date


def convert_ra_dec(ra, dec):
    # Convert RA from decimal degrees to HH MM SS.dd
    ra_angle = Angle(ra, unit=u.deg)
    ra_formatted = ra_angle.to_string(unit=u.hour, sep=' ', precision=2, pad=True)

    # Convert Dec from decimal degrees to sDD MM SS.d
    dec_angle = Angle(dec, unit=u.deg)
    dec_formatted = dec_angle.to_string(unit=u.deg, sep=' ', precision=1, alwayssign=True)

    return ra_formatted, dec_formatted



def format_to_MPC(matches, date, band, mag_list, file_path):
    """
    Reformats the given data into a file suitable for submission to the Minor Planet Center (MPC).
    https://www.minorplanetcenter.net/iau/info/ObsFormat.html

    Parameters:
    matches (list): A list of tuples containing matched objects from Pouakai and MPC.
    date (str): The observation date in the required format.
    band (str): The photometric band used for the observation.
    mag_list (list): List of observed magnitudes from the Pouakai side.
    savepath (str): The path where the file will be saved.

    Process:
    1. Extract the object designation from the MPC side of the matches.
    2. Extract the RA and DEC values from the Pouakai side of the matches.
    3. Extract the magnitude from the Pouakai side using mag_list.
    4. Format the extracted information into the required MPC format string.
    5. Save the formatted string to a file in the specified folder.

    Returns:
    None
    """
    observatory = '474'
    how_observed = 'C'   #MOA CCD image
    filename = file_path.split('/')[-1].split('.')[0]+"_MPC_submission.txt"
    formatted_date = convert_utc_format(str(date))
    full_path = os.path.join(os.getcwd(), "MPC_tests", filename)

    with open(full_path, 'w') as file:
        for match in matches:
            # Get object designation from matches (MPC side)
            object_designation = get_object_designation(matches)
        
            # Get RA and DEC from matches (Pouakai side)
            ra = match[0][3]
            dec = match[0][4]

            ra, dec = convert_ra_dec(ra, dec)  #Convert to the specified format

            # Get magnitude from mag_list (Pouakai side)
            mag = f'{mag_list[match[2]]:.1f}'
            

            # Format the data into the required MPC format string
            MPC_file = f"{object_designation}{'  '}{how_observed}{formatted_date}{' '}{ra:<12}{dec:<12}{' '*9}{mag}{' '}{band}{' '*6}{observatory:<3}"

            #To display the processes running
            print(MPC_file)

            # Write the formatted string to the file
            file.write(MPC_file + '\n')

    print(f"File saved to {full_path}")

    return


def MPC_maker(file_path):
    """
    Main function to process FITS files, extract header information, 
    and retrieve solar system objects for comparison and then 
    reformatting for the end result of a file suited for submission to the MPC.

    Steps:
    1. Extract header information from specified FITS file.
    2. Parse info to retrieve the list of solar system objects in the FoV.
    3. Compare this list with the source list from Pouakai.
    4. If compatibility is found, reformats all gathered info as a file for the MPC.

    Includes a try-except block to catch runtime errors.
    """
    try:
        info, data_pouakai = extracts_header_info(file_path)
        ra, dec, date, band, mag_limit, mag_list = info
        
        # Reformats the date into ISO format
        observation_time = '20' + date 

        # Get the list of Solar System objects in the FoV
        data_MPC = get_solar_system_objects(ra, dec, observation_time, mag_limit)

        # Compare this list with the source list from Pouakai
        tolerance = 0.1   #10 arcseconds is good!
        matches = compare_objects(data_pouakai, data_MPC, tolerance)

        filename = file_path.split('/')[-1].split('.')[0]+"_MPC_submission.txt"
        full_path = os.path.join(os.getcwd(), "MPC_tests", filename)

        if len(matches) !=0:
            print(f'{len(matches)} Match(es) Found!')
            # print(matches)
            format_to_MPC(matches, date, band, mag_list, file_path)
        else:
            with open(full_path, 'w') as file:
                file.write("None")

        print(f"File saved to {full_path}")

    except Exception as e:
        print(f"An error occurred: {e}")


# if __name__ == "__main__":
#     if len(sys.argv) != 2:
#         print("Usage: python process_fits.py <file_path>")
#     else:
#         file_path = sys.argv[1]
#         MPC_maker(file_path)


file_path1 = '/Users/redwy/MPC_Pouakai_pull_req/A6524-lsst2-R-6_phot.fits'
file_path2 = '/Users/redwy/MPC_Pouakai_pull_req/A6394-lsst1-R-10_phot.fits'
MPC_maker(file_path2)
