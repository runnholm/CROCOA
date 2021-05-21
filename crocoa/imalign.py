import numpy as np
import pandas as pd
import astrodrizzle as ad


def align(source_images, destination_dir, drizzle_config, temp_dir='./temp', reference_image=None, cleanup=True, hlet=False):
    pass


def match_images(
    reference_fits, image_fits, output_fits, verbose=True, normalization="white", clip=None
):
    """ Function that matches the coordinatesystem of image_fits
    to that of reference_fits by means of cross-correlation

    Parameters
    ----------
    reference_fits : str
        directory of the reference image
    image_fits : str
        directory of the image that should be matched
    output_fits : list
        path to original flc files where wcs should be updated
    """
    if verbose:
        print('ref', reference_fits)
        print('shift image', image_fits)

    # 1. Read data
    refdata = fits.getdata(reference_fits)
    refheader = fits.getheader(reference_fits)

    shiftdata = fits.getdata(image_fits)
    shiftheader = fits.getheader(image_fits)

    # 1.1 Normalize data
    if clip is not None:
        shiftdata = image_normalization(
            shiftdata, method="range", lower=0, higher=None,
        )
        refdata = image_normalization(
            refdata, method="range", lower=0, higher=None,
        )
    shiftdata_n = image_normalization(shiftdata, method=normalization)
    refdata_n = image_normalization(refdata, method=normalization)

    # 2. Get cross correlation shifts
    xshift, yshift = corr2d(shiftdata_n, refdata_n)
    if verbose:
        print("X shift: ", xshift)
        print("Y shift: ", yshift)

    # 3. Translate into Ra and dec
    xscale, yscale = get_pixel_scale(shiftheader)
    ra, dec = get_image_coords(shiftheader)
    dra, ddec = dpixels_to_dwcs(xshift, yshift, xscale, yscale, dec)
    if verbose:
        print("Delta RA: ", dra)
        print("Delta DEC: ", ddec)

    # 4. Write new coordinates to output
    for outfile in output_fits:
        if verbose:
            print("writing to: ", outfile)
        write_wcs_to_fits(outfile, dra, ddec)


def image_normalization(image, method="range", lower=0, higher=50):
    """ Normalizes an np.array

    Parameters
    ----------
    image : ndarray
        input image array
    method : {"white", "minmax", "range"}
        sets type of normalization. Default "white"

        minmax - subtract minimum value and divide by maximum

        white - subtract mean divide by standard deviation

        range - clips image data to range 0 to 100
    """
    if method is None:
        img = image
    elif method == "white":
        img = (image - image.mean()) / image.std()
    elif method == "minmax":
        img = (image - image.min()) / image.max()
    elif method == "range":
        img = np.clip(image, lower, higher)
    elif method == "clip_pixels":
        sorted_pixels = np.sort(image.flatten())
        lower_range = sorted_pixels[lower]

        higher_range = sorted_pixels[-higher]

        img = np.clip(image, lower_range, higher_range)
    elif method == "percentile":
        img = np.clip(image, np.percentile(image, lower),
                      np.percentile(image, higher))
    else:
        raise ValueError

    return img


def corr2d(image_to_shift, refimage):
    """Function for calculating the x and y shift required for image
    `image_to_shift` to be matched to `refimage`

    Parameters
    ----------
    image_to_shift : np.ndarray
        image_array
    refimage : np.ndarray
        image_array

    Returns
    -------
    x_shift : int
    y_shift : int
    """
    corr = correlate2d(image_to_shift, refimage, mode="same")
    shift = np.where(corr == corr.max())
    shift = (
        -1 * (shift[0] % refimage.shape[0] - refimage.shape[0] / 2 + 1),
        -1 * (shift[1] % refimage.shape[1] - refimage.shape[0] / 2 + 1),
    )

    return shift[1][0], shift[0][0]


def dpixels_to_dwcs(dx, dy, pixscale_x, pixscale_y, ref_dec):
    """ Converts from pixel shifts to shifts in ra and declination
    """
    d_dec = dy * pixscale_y

    d_ra = dx * pixscale_x * np.cos(np.deg2rad(ref_dec))
    return d_ra, d_dec


def get_pixel_scale(header):
    x_scale = header["CD1_1"]
    y_scale = header["CD2_2"]
    return x_scale, y_scale


def get_image_coords(header):
    """ Fetches RA and dec from a given header
    """
    return header["CRVAL1"], header["CRVAL2"]
