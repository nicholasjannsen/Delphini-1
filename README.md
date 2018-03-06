# Delphini-1
Software for the first AUSAT: Delphini-1. This software includes image reduction and aperture photometry code.

# Software descriptio. 
To explain the software in more detail we here present the results from it. The following code assumes that all data is placed  placed in the same folder, and the science, flat, bias, and dark frames are likewise called so (you can call them whatever you wnat you want, just remember to change the nmaes of the input files). The file 'test.py' can be used to make an easy test of the software. The following code example illustartes the usage:

```
from DELPHINI import DELPHINI
path = '/path/to/data/'
# Call class
XX = DELPHINI(path, 'science', 1, 1)

# Image Reduction:
XX.image_reduction('flat', 'bias')   

# Aperture photometry:
# Ellipse:
x_coor = [146, 201,  87, 213]
y_coor = [ 97,  52, 171, 208]
XX.aperture_photometry(x_coor, y_coor, ['ellipse', 6, 48, 8, 172], 'local')
# Trace:
x_coor = [107,  165,  49, 176]
y_coor = [103,   55, 176, 212]
XX.aperture_photometry(x_coor, y_coor, ['trace', 3, 78, 8, 172], 'local')
```

For the photometry software the first 2 entries in 'aperture_photometry' are the stellar coordinates. The next entry is the     "aperture" entry that takes 5 arguments: ['aperture', a, b, q, phi]. Here 'aperture' is either 'ellipse' or 'trace' corresponding to the two apertures. Because we are working with startrails, 'phi' is a tilt angle of the aperture between 0-180 degrees defined by the zero-point of the unit circle (hence counter clockwise from first quadrant). 

Using the ellipse as aperture "a" is the semi-minor axis of the ellipse, "b" is semi-major axis of the ellipse, "q" is the width of the local background flux from the ellipse. 

As mentioned the trace aperture is a mission-specific aperture that uses a circular aperture of radius "a". Given the coordinates of the most left part of the startrail, the COF is found inside this initial circular aperture. Next the circular aperture will be moved one pixel at a time in either the positive x or y direction depending on "phi"

```
x step: if  0<phi<45 or 135<phi<180  
y step: if 45<phi<135
```

For each pixel step in the x,y direction the opposite y,x pixel coordinate is determined by the COF. From our code example the aperture is moved in a x pixel step direction, which means for each step the belonging y coordinate is determine by the COF from the total circular aperture. Just as for elliptical aperture, "q" is here the width of the sky background aperture. The advantage of the trace aperture is, if it turns out that the satellite is very unstable, as long as the Signal to Noise Ratio (SNR) is sufficiently high, this routine will still follow the perhaps strange pattern of the COF for the stars.

The third argument for the utility "aperture_photometry" is if a local or global sky background flux should be used to correct the stellar flux. As mentioned above the local sky background flux is define by a band of width "q" around the stellar aperture. As the factor of stellar contamination and crowding is very hard to predict for our mission, the sky background flux can also be determined globally. This is done simply by slicing the image into s number of subframes. Inside each subframe n number of pixels having the lowest flux is found, hence, s*n is the total number of sky background pixels and the robust 3*median(sky-pixels)-2*mean(sky-pixels) value of the s*n number of pixels with lowest flux is then the sky background flux. When a high level of vignetting or other image artifacts is present the local sky background flux should be used. Also, the global sky background routine do not work for all subframes at the moment. 


# Output



# Usage and dependencies
Delphini-1 has been succesfully tested with python2.7 on Linux systems.

Some python packages need to be installed:

   1. numpy
   1. scipy
   1. pyplot
   1. matplotlib
   1. astropy
   1. PIL



All of the above can be installed using pip, e.g.:

```pip install numpy, scipy, matplotlib, pyfits, pycurl, ephem, rpy2```
