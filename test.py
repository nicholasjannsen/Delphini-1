"""
This is a utility to test the software in a easy and fast way. In the bottom is the main function, that takes
'name' as an input. Setting the name to e.g. 'image_reduction', the main function will call the 'test' 
function to calculate whatever that is in the function 'image_reduction'. First this an easy way to call 
utilities inside the class 'DELPHINI' and, secondly, it is easy to add extra small test by the 'if' statement.
To the this software from the terminal simply go to the directory of this file and type: './python test.py'.
"""
# Import functions:
import numpy as np

###########################################################################################################
#                                            TEST OF SOFTWARE                                             #
###########################################################################################################

# Call to class:
from Delphini_Photometry import Delphini_Photometry

# Several obtions can be tested separately:
def test(name):
     
     # Image Reduction:
     if name=='image_reduction':
          path = '/home/nicholas/Data/delphini-1/photometry/'
          XX = Delphini_Photometry(path, 'LF', 'fits', 1, 0)
          XX.image_reduction('FF', 'BF', 'DF')
          
     # Photometry:
     if name=='photometry':
          path = '/home/nicholas/Data/delphini-1/photometry/'
          XX = Delphini_Photometry(path, 'startrails', 'rgb', 1, 0)
          # Ellipse:
          x_coor = [146, 201,  87, 213]
          y_coor = [ 97,  52, 171, 208]
          XX.aperture_photometry(x_coor, y_coor, ['ellipse', 8, 48, 8, 172.5], 'local')
          # Trace:
          x_coor = [107,  165,  49, 176]
          y_coor = [103,   55, 176, 212]
          XX.aperture_photometry(x_coor, y_coor, ['trace', 3, 78, 8, 172.5], 'local')
     
if __name__ == '__main__': #---------------------------- Main function --------------------------#
     # CALL FUNCTIONS HERE BY CHANGIN NAME:
     test(name='photometry')    
