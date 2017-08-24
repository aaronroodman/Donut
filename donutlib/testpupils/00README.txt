README

The files in this folder provide 512x512 wavefront stamps from the DESI Zemax deck. All donuts are taken from the GFA sensor plane. Any file with "trim" in its name is readable by the python readimg function provided below. The file with the same name, but without "trim" includes information about the file from Zemax. For convenience, here are the file names and their corresponding x- and y- angles in degrees, in Zemax's coordinate system set up for the Echo 22 DESI deck.

fieldnames = ["center", "lower_left", "lower_right", "upper_left", "upper_right",
	      "field5", "field8", "field9", "field10", "field11", "field12"]

xangles = [1.318, #center
1.255, #lower left
1.286, #lower right
1.34,  #upper left
1.37,  #upper right
1.3,   #field5
1.34,  #field8
1.36,  #field9
1.3,   #field10
1.335, #field11
1.29   #field12]

yangles = [0.86,  #center
0.852,  #lower left
0.806,  #lower right
0.906,  #upper left
0.862,  #upper right
0.849,  #field5
0.86,   #field8
0.86,   #field9
0.873,  #field10
0.89,   #field11
0.83    #field12]

def readimg(filename, nbin):
    """Reads an example wavefront stamp in the "trim" format and turns it into a flat pupil image.

    Parameters
    ----------
    filename : string
        full path to the text file containing the wavefront stamp 
    nbin : int
        the number of pixels that the output stamp should have per side 


    Returns
    -------
    img : array
        the final pupil stamp with the number of pixels per side defined by nbin
    """
    import numpy as np
    import pandas as pd
    from skimage.transform import resize

    img = pd.read_csv(filename, sep='\s+', header=None)
    img = img.as_matrix()
    img[img != 0] = 1
    img = resize(img, (nbin, nbin), mode='reflect')
    return img
