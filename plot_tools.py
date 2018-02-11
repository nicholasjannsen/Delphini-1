""" 
This utility was made to hide away routines that plot the data and thus making the code in the DELPHINI class
looks more simple. Also for future work more routines can be added into this utility. For now it only have a 
function 'FITS' that easily plot fits or RGB images, and a function that plot ellipses 'plot_ellipse'.
"""

# Import numpy:
import numpy as np
import matplotlib.pyplot as plt
import math

###########################################################################################################
#                                            UTILITY FUNCTIONS                                            #
###########################################################################################################

def FITS(img, scale, ds=2):
    """
    This function makes it easier to plot fits/jpg/png etc. images. Its uses the 'Image_Scale' utility to set
    the scale of the background level. To plot something on top of this image, one needs to apply 'plt.show()'
    after 'FITS()' is called.
    ----------INPUT:
    img            : Image that you want to plot
    scale          : Scale of background level
    ds             : This how many standard deviations the level should be determined with.
    """
    import pylab
    from image_scale import linear, sqrt, log, asinh
    if scale=='linear': scale=linear
    if scale=='sqrt'  : scale=sqrt
    if scale=='log'   : scale=log
    if scale=='asinh' : scale=asinh
    img_min, img_max = img.mean()-ds*img.std(), img.mean()+ds*img.std()
    pylab.clf(); pylab.imshow(scale(img, img_min, img_max), cmap='gray', origin='lower')


def plot_ellipse(semimaj=1, semimin=1, phi=0, x_cent=0, y_cent=0, mark='b-'):
    """
    An easy to use function for plotting ellipses in Python 2.7! The function creates a 2D ellipse in 
    polar coordinates then transforms to cartesian coordinates. It can take a covariance matrix and plot 
    contours from it.
    ----------INPUT:
    semimaj        : length of semimajor axis
    semimin        : length of semiminor axis
    phi            : angle in radians of semimajor axis above positive x axis
    x_cent         : X coordinate center
    y_cent         : Y coordinate center
    """
    # Convert degress to radias:
    phi = math.radians(phi)
    # Generate data for ellipse structure
    theta = np.linspace(0,2*np.pi,1e3)
    r = 1 / np.sqrt((np.cos(theta))**2 + (np.sin(theta))**2)
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    data = np.array([x,y])
    S = np.array([[semimaj,0],[0,semimin]])
    R = np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]])
    T = np.dot(R,S)
    data = np.dot(T,data)
    data[0] += x_cent
    data[1] += y_cent
    # Plot:
    plt.plot(data[0], data[1], mark)
