####### Import the necessary modules ############
import sim_SAM_SPHERE
import pdb
import calibrate_SAM
import os


#######  The tolowing parameters are necessary to run the SAM pipelie ##################
mask_filename = '7holes_jwst_mask.txt' #A text file with the mask geometry
wave = 3.828e-06 #4.2934e-6 #Central wavelength of the observations (meters)
bandwidth = 0.205e-06 #0.202e-6 #Bandwidth of the observations (meters)
hole_size = 0.82
imsize = 80
px_scale = 65.6
hole_geom = 'HEXAGON'
inst = 'JWST'
arrname = 'SIM' ## This could be DATA or SIM
rotation_angle = 0.0 ### In case the msk is not propely aligned with the position indicated in the manual
oversample = 1.0 ## Number of times that you want to oversample the data

data_filename = 'sci_psf3_w1.txt'
source = 'disk'
sim_SAM_SPHERE.simSAM_PSF (data_filename, mask_filename,  wave, bandwidth, hole_size, px_scale, imsize, hole_geom, source, inst, \
           arrname, rotation_angle, oversample)


#######  The tolowing parameters are necessary to run the SAM pipelie ##################
mask_filename = '7holes_jwst_mask.txt' #A text file with the mask geometry
wave = 4.286e-06 #4.2934e-6 #Central wavelength of the observations (meters)
bandwidth = 0.202e-06 #0.202e-6 #Bandwidth of the observations (meters)
hole_size = 0.82
imsize = 80
px_scale = 65.6
hole_geom = 'HEXAGON'
inst = 'JWST'
arrname = 'SIM' ## This could be DATA or SIM
rotation_angle = 0.0 ### In case the msk is not propely aligned with the position indicated in the manual
oversample = 1.0 ## Number of times that you want to oversample the data

data_filename = 'sci_psf3_w2.txt'
source = 'disk'
sim_SAM_SPHERE.simSAM_PSF (data_filename, mask_filename,  wave, bandwidth, hole_size, px_scale, imsize, hole_geom, source, inst, \
           arrname, rotation_angle, oversample)


#######  The tolowing parameters are necessary to run the SAM pipelie ##################
mask_filename = '7holes_jwst_mask.txt' #A text file with the mask geometry
wave = 4.817e-06 #4.2934e-6 #Central wavelength of the observations (meters)
bandwidth = 0.298e-06 #0.202e-6 #Bandwidth of the observations (meters)
hole_size = 0.82
imsize = 80
px_scale = 65.6
hole_geom = 'HEXAGON'
inst = 'JWST'
arrname = 'SIM' ## This could be DATA or SIM
rotation_angle = 0.0 ### In case the msk is not propely aligned with the position indicated in the manual
oversample = 1.0 ## Number of times that you want to oversample the data

data_filename = 'sci_psf3_w3.txt'
source = 'disk'
sim_SAM_SPHERE.simSAM_PSF (data_filename, mask_filename,  wave, bandwidth, hole_size, px_scale, imsize, hole_geom, source, inst, \
           arrname, rotation_angle, oversample)
