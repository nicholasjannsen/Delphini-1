"""
CLASS DESCRIBTION:
------------------
This is the astronomical pipline used to reduce and obtain photometry for the AUSAT (the first satellite 
lunched by Aarhus University). This pipeline was constructed as a final project for the 2nd workshop of 
Delphini-1. The overall aim for the routines in this DELPHINI python class, is to perform 'image reduction' 
and 'photmetry'.To each routine a description for its functionality is given. For the functionality see the
README file downloaded alongside with this class To get this class working the utilities 'plot_tools.py' and 
'image_scale.py' needs to be placed in the same folder (or you can include a
path to the utilities and it will work as well).
"""

# Numpy: 
import numpy as np
from numpy import sum
from numpy import inf, nan, sin, cos, tan, arctan, pi, sqrt, diff, std, diag, argmin, log10, meshgrid
from numpy import mean, median, nanargmax, zeros, ones, ceil, delete, shape, roll, nonzero, std, power
from numpy import arange, array, asarray, size, vstack, hstack, copy, loadtxt, where, savetxt, linspace
from numpy import shape, copy, sort, histogram, meshgrid, ogrid, arctan2, dot
# Packages:
import time, sys, glob 
import pyfits, pylab, scipy, math
# Functions:
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter
import scipy.ndimage as snd
from scipy.misc import imsave
from astropy.io import fits
from astropy.time import Time
from PIL import Image
# Own functions:
from plot_tools import FITS

###########################################################################################################
#                                            DEFINE CLASS                                                 #
###########################################################################################################

class Delphini_Photometry(object):
    # INITILIZE THE CLASSE: 
    def __init__(self, path, LF_name, imgtype, plot, save):
        
        # DEFINE GLOBAL VARIABLES (DGV)

        # Customized informations:
        self.path    = path       # Directory path to data
        self.LF_name = LF_name    # Name of Light Frame (LF)
        self.imgtype = imgtype    # Format of images
        self.plot    = plot       # Plot if plot==1
        self.save    = save       # Save if save==1

        # Load light frames:
        self.LF_files = sort(glob.glob('{}{}*'.format(self.path, self.LF_name)))
        if self.imgtype=='fits' :
            self.LF_i = array([pyfits.getdata(str(files)) for files in self.LF_files])
            hdulist = fits.open('{}'.format(str(self.LF_files[0])))
            self.t_exp_LF = hdulist[0].header['EXPTIME']            # Exposure time (CHANGE HEADER NAME!!!)
        if self.imgtype=='rgb': # Here we combine the colors with 'sum, axis=3':
            self.LF_i = sum(array([asarray(Image.open(str(files))) for files in self.LF_files]), axis=3)
        
        # Image dimensions:
        self.n, self.h, self.w = shape(self.LF_i)

        # HEADER INFORMATION:
        # load times:
        
        #time[i] = Time(hdulist[0].header['date'], scale='utc').jd  # Scaling to utc time and Julian Date
        #self.DF = 1 # Dark value from header:    
        
###########################################################################################################
#                                            CLASS FUNCTIONS                                              #
###########################################################################################################
        
    def image_reduction(self, FF_name, BF_name, DF_name):
        """ 
        This routine is a tool to reduce and correct for image artifacts such as noise in the optics, 
        heat noise, dust on the lens and mirrors, vignetting etc.. This function loads all available Light-
        Frames (LF), Flat-Frames (FF), and Bias-Frames (BF) from a folder given by directory path. For our 
        satellite the flats and bias' are obtained before lunch, however, the Dark-Frame (DF) median value 
        which is given in the image header. 
        ----------INPUT:
        path           : Directory path to data.
        LF_name (self) : Name of Light Frames (LF) except end number.
        DF_value (self): Median value from dark frame.
        FF_name        : Name of Flat Frames (FF) except end number.
        BF_name        : Name of Bias Frames (BF) except end number.
        DF_name        : Name of Dark Frames (DF) except end number.
        plot           : If plot=1: plots all relevant frames.
        save           : If save=1: saves all corrected light frames in seperate files.
        ---------OUTPUT:
        CF             : Corrected Frames (CF). Can be also be saved. """
        print '--------------------------------------------------------------------- image_reduction'
        start_time = time.time()  # Take time

        # Construct 3D cube of light-frames and flat-frames with dim(number, x-dim, y-dim):
        FF_files = sort(glob.glob('{}{}*'.format(self.path, FF_name)))
        BF_files = sort(glob.glob('{}{}*'.format(self.path, BF_name)))
        DF_files = sort(glob.glob('{}{}*'.format(self.path, DF_name)))

        # Load images into an array:
        if self.imgtype=='fits':
            FF_i = array([pyfits.getdata(str(files)) for files in FF_files])
            BF_i = array([pyfits.getdata(str(files)) for files in BF_files])
            DF_i = array([pyfits.getdata(str(files)) for files in DF_files])
        if self.imgtype=='rgb' : # Combine RGB:
            FF_i = sum(array([asarray(Image.open(str(files))) for files in FF_files]), axis=3)
            BF_i = sum(array([asarray(Image.open(str(files))) for files in BF_files]), axis=3)
            DF_i = sum(array([asarray(Image.open(str(files))) for files in DF_files]), axis=3)  
            
        # Median correction images:
        FF = median(FF_i, axis=0) 
        BF = median(BF_i, axis=0)
        DF = median(DF_i, axis=0)
        
        # Median correction images:
        BF = median(BF_i, axis=0)                # Master bias 
        # Correct for bias:
        FF_m_b = median(FF_i-BF, axis=0)
        FF_i_b = FF_i - BF
        DF_b_m = median(DF_i-BF, axis=0)         # Dark Current
        # Scale dark to science frames:
        hdulist  = fits.open('{}'.format(str(DF_files[0])))
        t_exp_DF = hdulist[0].header['EXPTIME']  # Exposure time (CHANGE HEADER NAME!!!)
        DF = self.t_exp_LF/t_exp_DF * DF_b_m     # Master dark   
        # Correct flats:
        FF = median(FF_i_b - DF, axis=0)         # Master flat
        # Final correction:
        CF_i = (self.LF_i - BF - DF)/FF          # Calibrated science images       
        print 'Filter done in time: %0.5f s' %(time.time()-start_time)
        
        # Plot if you like:
        if self.plot==1:
            FITS(FF, 'linear', 3); plt.show() 
            FITS(BF, 'linear', 3); plt.show()
            FITS(DF, 'linear', 3); plt.show()   
            FITS(self.LF_i[0], 'linear', 2); plt.show()
            FITS(CF_i[0], 'linear', 2); plt.show()
    
        # Save if you like:
        if self.save==1:
            [imsave('{}CF_%02d.{}'.format(self.path, self.imgtype) %i, CF_i[i]) for i in range(self.n)]
            

            
    def aperture_photometry(self, x_stars, y_stars, aperture, background):
        """ 
        This routines perform aperture photometry for stellar images with stars as trails. This function
        returns flux for all the stars that asigned coordinates as well as the Signal to Noise Ratio (SNR)
        for these stars. 
        ----------INPUT:
        x_star         : x stellar coordinate
        y_star         : y stellar coordinate
        aperture       : Aperture to be used: ['ellipse' or 'trace', a, b, q, phi]
        background     : Sky background to be used; 'local' or 'global' 
        """
        print '--------------------------------------------------------------------- aperture_photometry'

        #--- CONSTANTS ---#
        gain = 0.73          # Gain of camera: electrons pr ADU (ADU = counts from object)- ajust to camera!
        ron  = 3.3           # Read out noise - ajust to camera!
        con  = 25            # Magnitude constant

        #--- PHOTOMETRY ---#
        # Find fluxes:
        N = len(x_stars) # Number of stellar objects 
        flux_star = zeros((self.n,N))
        SNR_i = zeros((self.n,N))
        for i in range(self.n):  # Loop over all images: if timeseries is available
            for j in range(N):   # Loop over all stars and find flux: using same aperture size
                flux_sky, n_star_pix, flux_star[i][j] = self.aperture(self.LF_i[i], x_stars[j], y_stars[j],\
                                                                      aperture, background)
                SNR_i[i][j] = self.SNR(flux_sky, n_star_pix, flux_star[i][j], gain, ron)

        #--- FINAL CORRECTIONS ---#
        print flux_star, flux_sky, SNR_i
        return flux_star, SNR_i
            

                
    def aperture(self, LF, x, y, aperture, background):
        """
        This function calculate the stellar and sky background flux either using a highly elliptical aperture
        or a aperture that traces the the Center Of Flux (COF) to fit the startrails. The routine returns the
        stellar flux 'flux_star', the sky background flux 'flux_sky', and the number of stellar pixels
        'n_pix_star'.
        ----------INPUT:
        LF             : A single Light Frame (LF)
        x, y           : Stellar coordinate
        aperture       : aperture is either: ['ellipse' or 'trace', a, b, q, phi]
        background     : sky background is either: 'local' or 'global'
        """
        # Aperture only handle integers:
        x, y = int(x), int(y)
        
        # Grid for aperture:
        x_grid, y_grid = meshgrid(range(self.w), range(self.h))
        grid = LF[y_grid, x_grid]  # Grid
        xx = x_grid - x            # Displacement for stellar x coor
        yy = y_grid - y            # Displacement for stellar y coor
        
        # Ellipse parameters:
        a   = aperture[1]                # Semiminor/lenght axis
        b   = aperture[2]                # Semimajor/width axis 
        q   = aperture[3]                # Width of background band
        phi = math.radians(aperture[4])  # Tilt angle: [0:180] deg
        
        # Axis of background: 
        a_sky = a + q   # Semiminor/lenght axis
        b_sky = b + q   # Semimajor/width axis

        # Global background:
        if background=='global':
            flux_sky = self.global_sky_background(LF)
        
        #--- ELLIPTIC APERTURE ---#
        if aperture[0]=='ellipse':

            # Parametrasation and rotation of the ellipse (1 = x^2/a^2 + y^2/b^2)
            phi = phi - pi/2 # The ellipse have an offset from unity circle 
            EE_star  = (xx*cos(phi)+yy*sin(phi))**2/a**2     + (xx*sin(phi)-yy*cos(phi))**2/b**2
            EE_sky   = (xx*cos(phi)+yy*sin(phi))**2/a_sky**2 + (xx*sin(phi)-yy*cos(phi))**2/b_sky**2
            
            # Local background:
            if background=='local':
                E_sky    = ((EE_star>1)*(EE_sky<1))*grid   # Sky background determined by width q 
                sky      = E_sky[nonzero(E_sky)]           # Sky pixel values
                flux_sky = 3*median(sky) - 2*mean(sky)     # Robust sky background flux
            
            # Star: 
            E_star     = (EE_star<=1)*grid
            star       = E_star[nonzero(E_star)]-flux_sky  # Stellar corrected pixels
            n_pix_star = sum(EE_star<=1)                   # Number of used star pixels
            flux_star  = sum(star)                         # Flux from star
            
        #--- TRACE APERTURE ---#
        if aperture[0]=='trace':

            # Star initial:
            CC_star  = sqrt(xx**2 + yy**2) - a      # Start frame used to trace from:
            star_img = (CC_star<=1)*grid            # Star images
            star_val = star_img[nonzero(star_img)]  # Stellar pixels
            n_pix_star = len(star_val)              # Number of pixels in circle
                
            # Loop-step in x or y depends on phi:
            # Step in x and finds y centroid for each step:
            if 0<=phi<=pi/4 or pi*3/4<=phi<=pi:
                step = 'x step'
                x_step = range(b)
                y_step = zeros(b)
            # Step in y and finds x centroid for each step:
            if pi/4<phi<pi*3/4:
                step = 'y step'
                x_step = zeros(b)
                y_step = range(b)

            # Loop over trace step:
            x_cen = zeros(b)
            y_cen = zeros(b)
            star_true = zeros((b, self.h, self.w))
            sky_true  = zeros((b, self.h, self.w))
            for i in range(b):
                # Find Center Of Flux (COF):
                y_cen[i], x_cen[i] = self.center_of_flux(star_img, n_pix_star)
                if step=='x step': XX = x + x_step[i]; YY = y_cen[i] + y_step[i]
                if step=='y step': YY = y + y_step[i]; XX = x_cen[i] + x_step[i]
                # A circular aperture is used to trace with:
                CC_star = sqrt((x_grid-XX)**2 + (y_grid-YY)**2) - a
                star_true[i] = (CC_star<=1)                      # Array of true and false statement: star
                star_img     = star_true[i]*grid                 # Image: star
                n_pix_star   = len(star_img[nonzero(star_img)])  # Number of pixels: star
                # If local background:
                if background=='local':
                    CC_sky  = sqrt((x_grid-XX)**2 + (y_grid-YY)**2) - a_sky
                    sky_true[i] = ((CC_star>1)*(CC_sky<1))*grid  # Array of true and false statement: sky

            # Local sky flux:
            if background=='local':
                star_false = np.logical_not(np.sum(star_true, axis=0))  # Invert to get rid of star
                sky_true_x = np.sum(sky_true, axis=0) > 0
                sky_img  = sky_true_x.astype(np.int)*star_false*LF      # Here bol: True*False = False
                sky      = sky_img[nonzero(sky_img)]
                flux_sky = 3*median(sky) - 2*mean(sky)                  # Robust sky flux
            
            # Stellar flux: 
            star_true_x = np.sum(star_true, axis=0) > 0         # Bool array of stellar pixels
            star_img    = star_true_x.astype(np.int)*LF         # Stellar image
            star        = star_img[nonzero(star_img)]-flux_sky  # Stellar corrected pixels
            n_star_pix  = len(star)
            flux_star   = sum(star_img)
            
        # PLOT IF YOU LIKE::
        if self.plot==1:
            # Elliptic aperture with local background:
            if aperture[0]=='ellipse':
                from plot_tools import plot_ellipse
                FITS(LF, 'linear', 2)
                plot_ellipse(a,     b,     math.degrees(phi), x, y, 'g')  # Stellar aperture
                plot_ellipse(a_sky, b_sky, math.degrees(phi), x, y, 'm')  # Background aperture
                plt.show()
            # Box aperture with local background:
            if aperture[0]=='trace':
                 t = linspace(0, 2*pi)
                 FITS(LF, 'linear', 2)
                 [plt.plot(a_sky*cos(t)+(x+x_step[i]), a_sky*sin(t)+y_cen[i], 'b-') for i in range(b)]
                 [plt.plot(a*cos(t)+(x+x_step[i]), a*sin(t)+y_cen[i], 'g-') for i in range(b)]
                 plt.show()

        return flux_sky, n_pix_star, flux_star

    

    def global_sky_background(self, LF):
        """
        This function calculates the global sky background flux. 200 minimum pixels are used as background
        as a standard. As the background flux may not be uniform across CCD in space, the light frame is
        divided into 's' number of subframes, thus, in each subframe a total number of 'n' number of pixels
        are found across the sky and used in the sky background flux.
        ----------INPUT:
        LF             : A single Light Frame (LF)
        ---------OUTPUT:
        flux_sky       : Global sky background flux
        """
        # Variables:
        s     = 9                             # Number of subframes (CHANGE IF NEEDED!) E.g. 4, 9, 16 etc. 
        n     = self.h*self.w/(self.h+self.w) # Number of pixels used in subframes scales with image dim 
        nrows = self.h/(s/2)                  # Numbers of rows in each subframe
        ncols = self.w/(s/2)                  # Numbers of columns in each subframe

        # Reshape light frame into subframe:
        LF_sub = (LF.reshape(self.h//nrows, nrows, -1, ncols).swapaxes(1,2).reshape(-1, nrows, ncols))

        # Loop over all subframes:
        min_val = np.zeros((s,n))
        for i in range(s):
            # Loops over all pixels:
            for j in range(n):
                min_val[i,j] = np.min(LF_sub[i])                # Minimum value for array
                min_dex = np.where(LF_sub[i] == min_val[i,j])   # Find row, column for min value
                # Min pixel is set to max in order to find the next min:
                LF_sub[i, min_dex[0][0], min_dex[1][0]] = np.max(LF_sub[i]) 

        # Flux:
        flux_sky = 3*median(min_val) - 2*mean(min_val)  # Mean flux from pixels
        return flux_sky


    
    def SNR(self, flux_sky, n_pix_star, flux_star, gain, ron):
        """
        This function calculates the Signal-to-Noise Ratio (SNR). 
        ----------INPUT:
        flux_sky       : Flux from each sky pixel
        n_pix_star     : Number of pixels used to find star flux
        flux_star      : Flux from stellar object
        gain           : Gain CCD (e-/ADU)
        ron            : Read out noise (ron) (e-)
        """
        SNR = (gain*flux_star/sqrt(gain*flux_star + n_pix_star*gain*flux_sky + n_pix_star*ron**2))      
        return SNR


    
    def center_of_flux(self, LF, n_pix):
        """
        This function finds the center of flux for all desired stellar object. Here LF is the masked image
        thus every pixel is set to zero except the star) and n_pix is the number of pixels one wish to use
        in order to find the COF.
        """
        # Loops over all pixels:
        LF_copy  = copy(LF)     # Copy to avoid overwriting
        flux_max = zeros(n_pix)
        x_max = zeros(n_pix)
        y_max = zeros(n_pix)
        pixel = zeros(n_pix)
        for j in range(n_pix):
            flux_max[j] = np.max(LF_copy)               # Maximum value for array
            max_dex = np.where(LF_copy == flux_max[j])  # Find row, column for min value
            x_max[j] = max_dex[0][0]                    # max for x coordinate
            y_max[j] = max_dex[1][0]                    # max for y coordinate
            pixel[j] = j
            # Min pixel is et to max in order to find the next min:
            LF_copy[int(y_max[j]), int(x_max[j])] = 0 

        # Flux center is found:
        flux = sum(pixel)
        cen_x = 1/flux*dot(pixel,x_max)              # Flux-center in x
        cen_y = 1/flux*dot(pixel,y_max)              # Flux-center in y
        return cen_x, cen_y
