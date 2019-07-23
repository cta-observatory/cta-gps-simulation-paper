import numpy as np
from astropy.io import fits


def normalize_template(in_fits, in_hdu, out_fits,
                       thresh=1.e-6):
    """
    Normalize a 2D WCS map so that integral over solid angle = 1.
    Implementation only works for small distortions, i.e.,
    for small map sizes and TAN projection
    :param in_fits: name of input fits file
    :param in_hdu: input hdu name or number
    :param out_fits: name of output fits file
    :param thresh: threshold of map value for output map as fraction of max value
    :return:
    """

    # input hdu
    in_hdu = fits.open(in_fits)[in_hdu]

    # test if the input map has 2 dimensions
    if len(np.shape(in_hdu.data)) == 2:
        pass
    else:
        print("ERROR: template normalization only implemented for 2D maps")

    # test if projection distortion can be ignored, here I use 3 deg as max size
    if 'TAN' in in_hdu.header['CTYPE1']\
            and np.abs(in_hdu.header['NAXIS1'] * in_hdu.header['CDELT1']) < 3.\
            and np.abs(in_hdu.header['NAXIS2'] * in_hdu.header['CDELT2']) < 3.:
        pass
    else:
        print("ERROR: template normalization only works for small maps in TAN projection")

    # filter low-value pixels
    in_hdu.data[in_hdu.data < thresh * np.max(in_hdu.data)] = 0.

    # normalize map
    in_hdu.data /= np.sum(in_hdu.data)

    # write normalized map to disk
    in_hdu.writeto(out_fits)

    return
