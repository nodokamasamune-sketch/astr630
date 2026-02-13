############################################################################
#
# plot_xsl_dr3.py
#
#
# Description:
#       Read the binary DR3 XSL spectra
#
# Functions:
#       -- read_bin_spec: Read a binary spectrum
#       -- example
#
#
# History:
#       - 27 / 11 / 19 : Creation of plot_xsl_dr2.py by A. Gonneau
#       - 17 / 02 / 22 : Modifying for DR3 data plot_xsl_dr3.py by K.Verro
#
#
############################################################################

## Load some Python modules
from astropy.io import fits
import numpy as np


############################################################################

# **********************************************************
#
# 	Read a binary spectrum
#
# **********************************************************

def read_dr3_spec(spec_name, ang='', star_kw='', loss_corr_kw='',  av_kw=''):


	# If the file exists
    try: 
	
        ###########

        # Print some info
        print(' ')
        print('File to read:', spec_name)


		###########

		# Open the spectrum
        hdu = fits.open(spec_name)
		

		###########

		# Define the columns based on file name
        keyword1 = '_ncl.'
        keyword2 = '_ncge.'
        keyword3 = '_scl.'
       
        if keyword1 in spec_name:
            print('This spectrum is NOT corrected for slit flux losses.') 
            flux = hdu[1].data['FLUX']
            flux_dr = hdu[1].data['FLUX_DR']
            waves = hdu[1].data['WAVE']
            flux_dr_final = np.array(flux_dr)
            
        elif keyword2 in spec_name:
            print('spectrum IS NOT corrected for galactic dust extinction.')
            print('De-reddened spectra unavailable.')
            flux = hdu[1].data['FLUX']
            waves = hdu[1].data['WAVE']
            flux_dr_final = np.empty(len(waves))
            flux_dr_final[:] = np.NAN
        elif keyword3 in spec_name:
            print('spectrum is corrected for galactic dust extinction' )
            print('and flux losses with a spline function')
            flux = hdu[1].data['FLUX']
            flux_dr = hdu[1].data['FLUX_SC']
            waves = hdu[1].data['WAVE']
            flux_dr_final = np.array(flux_dr)
        else:
            flux = hdu[1].data['FLUX']
            flux_dr = hdu[1].data['FLUX_DR']
            waves = hdu[1].data['WAVE']
            flux_dr_final = np.array(flux_dr)	
            
		###########

        # Convert the arrays to numpy arrays 
        waves_final = np.array(waves)
        flux_final = np.array(flux)
        	

		###########

		# Get the units of the waves
        prim_hdr = hdu[0].header
        sec_hdr = hdu[1].header

        unit = ''
        unit = sec_hdr['TUNIT1']
		

		###########

		# Change the units to Angstrom
        if (len(ang) != 0):
            if (unit == 'nm'):
                unit = 'A'
                waves_final = waves_final * 10.



		###########

        # Get the star name
        if (len(star_kw) != 0):
            star = ''
            star = prim_hdr['HNAME']
            
        ###########

        # Get the Av of dust extinction correction
        # and the origin of that value
        if (len(star_kw) != 0):
            av = ''
            av = prim_hdr['AV_VAL']
            av_ori = ''
            av_ori = prim_hdr['AV_ORI'] 
            
		###########

        # Get the flux-loss kw
        if (len(loss_corr_kw) != 0):
            loss_corr = ''
            loss_corr = prim_hdr['LOSS_COR']
            if loss_corr == False:
                loss_corr_spline = prim_hdr['SPL_COR']
                if loss_corr_spline == True:
                    loss_corr = 'spline'
            

	###########

        # File not found 
    except IOError:
        print('=> File not found: ', spec_name)
        flux_final = 999
        waves_final = 999
        unit = '999'


	###########

        # Define the final results
    list_final = [flux_final, flux_dr_final, waves_final, unit]
  
	# Return the flux, waves, unit
    if (len(star_kw) != 0):
        list_final.append(star)

    if (len(loss_corr_kw) != 0):
        list_final.append(loss_corr)
    
    if (len(av_kw) != 0):
        list_final.append(av)
        list_final.append(av_ori)


    # Return the results
    return list_final
	
#example:	
# spec_name = "xsl_spectrum_X0703_merged.fits"
# [flux_temp, flux_tempdr_temp, waves_temp, unit, star, loss_corr, Av, Av_ori]= read_dr3_spec(spec_name, ang='on', star_kw='on', loss_corr_kw='on', av_kw = 'on')
# print([flux_temp, flux_tempdr_temp, waves_temp, unit, star, loss_corr, Av, Av_ori])
