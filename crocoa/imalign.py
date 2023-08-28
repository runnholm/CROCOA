import numpy as np
import pandas as pd
import glob
import os
from drizzlepac.astrodrizzle import AstroDrizzle as adriz
from crocoa import filemanagement as fm
from crocoa.utilities import suppress_stdout_stderr
from astropy.io import fits
from scipy.signal import correlate2d
from skimage.registration import phase_cross_correlation
from datetime import datetime

def match_images(
    reference_fits,
    image_fits,
    verbose=True,
    normalization="white",
    clip=None,
    mask_method="zero",
):
    """Function that matches the coordinatesystem of image_fits
    to that of reference_fits by means of cross-correlation

    Parameters
    ----------
    reference_fits : str
        directory of the reference image
    image_fits : str
        directory of the image that should be matched
    """
    if verbose:
        print("ref", reference_fits)
        print("shift image", image_fits)

    # 1. Read data
    refdata = fits.getdata(reference_fits)
    refheader = fits.getheader(reference_fits)

    shiftdata = fits.getdata(image_fits)
    shiftheader = fits.getheader(image_fits)

    # 1.01 Check if there are nans in data and create masks in that case
    # masks are True on valid pixels
    if np.any(np.isnan(refdata)):
        refmask = ~np.isnan(refdata)
    else:
        refmask = np.ones(refdata.shape, dtype="bool")
    if np.any(np.isnan(shiftdata)):
        shiftmask = ~np.isnan(shiftdata)
    else:
        shiftmask = np.ones(shiftdata.shape, dtype="bool")

    # 1.1 Normalize data
    if clip is not None:
        shiftdata = image_normalization(
            shiftdata,
            method="range",
            lower=0,
            higher=None,
            mask=shiftmask,
        )
        refdata = image_normalization(
            refdata,
            method="range",
            lower=0,
            higher=None,
            mask=refmask,
        )
    shiftdata_n = image_normalization(shiftdata, method=normalization, mask=shiftmask)
    refdata_n = image_normalization(refdata, method=normalization, mask=refmask)

    # 2. Get cross correlation shifts
    if (np.any(np.isnan(refdata)) | np.any(np.isnan(shiftdata))) & (
        mask_method == "skimage"
    ):
        # run scikit correlate to deal with masked pixels
        print("Warning: NaN pixels in input data (skimage masked correlation chosen).")
        xshift, yshift = corr2d_sk(refdata_n, shiftdata_n, refmask, shiftmask)
    else:
        shiftdata_n = np.where(shiftmask, shiftdata_n, 0.0)
        refdata_n = np.where(refmask, refdata_n, 0.0)
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

    return dra, ddec


def align_multiple_filters(image_sets, reference_set_index=0, cleanup=True, matching_config={}, perform_manual_shifts=True):
    """ Function for aligning image sets between multiple filters 
    Parameters
    ----------
    image_sets : list or iterable
        list of ImageSet instances
    reference_set_index : int
        index of the ImageSet in the image_set list that is to be used as reference for cross correlation
    cleanup : bool
        whether to cleanup the temporary folders and intermediate drizzling results
    matching_config : dict
        keyword arguments that are passed directly to match images
    perform_manual_shifts : bool
        whether or not to apply a manual shift. This manual shift should be given to the ImageSet 
        during construction
    """
    source_images = []
    for image_set in image_sets:
        image_set.make_all_copies()
        if perform_manual_shifts:
            image_set.apply_manual_shifts()
        image_set.drizzle()
        source_images.append(image_set.drizzled_files)
    
    reference_image = source_images[reference_set_index]
    for i, image in source_images:
        if i == reference_set_index:
            pass
        else:
            dra, ddec = match_images(reference_image, image, **matching_config)
            image_sets[i].backpropagate_wcs_shift(dra, ddec)
    if cleanup:
        image_set.clean_temp_directories()

def align_single_filter(image_set, reference_image_index=0, cleanup=True, matching_config={}, perform_manual_shifts=True):
    """ Function for aligning images in a single set
    Parameters
    ----------
    image_set : ImageSet instance
        an ImageSet instance with the relevant files to be aligned
    reference_image_index : int
        index of the image in the image_set list that is to be used as reference for cross correlation
    cleanup : bool
        whether to cleanup the temporary folders and intermediate drizzling results
    matching_config : dict
        keyword arguments that are passed directly to match images
    perform_manual_shifts : bool
        whether or not to apply a manual shift. This manual shift should be given to the ImageSet 
        during construction
    """
    image_set.make_all_copies()
    if perform_manual_shifts:
        image_set.apply_manual_shifts()
    image_set.drizzle(individual=True)
    source_images = image_set.drizzled_files
    reference_image = source_images[reference_image_index]
    for i, image in enumerate(source_images):
        if i == 0:
            pass
        else:
            dra, ddec = match_images(reference_image, image, **matching_config)
            image_set.images[i].backpropagate_wcs_shift(dra, ddec)
    
    if cleanup:
        image_set.clean_temp_directories()




def manual_wcs_shift(image_list, shift_dict):
    """Function for applying a manual alteration of the WCS
    Primarily for doing rough alignment before correlation analysis
    """
    for img in image_list:
        visit_name = os.path.basename(img).split("_")[0]
        if visit_name in shift_dict.keys():
            fm.write_wcs_to_fits(
                img, shift_dict[visit_name]["dra"], shift_dict[visit_name]["ddec"]
            )
            print("Applied manual correction to: ", img)


def image_normalization(image, method="range", lower=0, higher=50, mask=None):
    """Normalizes an np.array

    Parameters
    ----------
    image : ndarray
        input image array
    method : {"white", "minmax", "range"}
        sets type of normalization. Default "white"

        minmax - subtract minimum value and divide by maximum

        white - subtract mean divide by standard deviation

        range - clips image data to range 0 to 100
    mask : ndarray (boolean)
        mask array for input image, valid pixels are True.

        Default: All pixels valid. NB masked pixels are set to 0 for statistics.
    """

    if mask is None:
        mask = np.ones(image.shape, dtype="bool")

    image = np.where(mask, image, 0.0)
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
        img = np.clip(image, np.percentile(image, lower), np.percentile(image, higher))
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


def corr2d_sk(refima, ima, mask_ref, mask_ima):
    """Performs (masked) cross correlation with scikit_image"""
    shift = phase_cross_correlation(
        refima, ima, reference_mask=mask_ref, moving_mask=mask_ima
    )
    # shifts are given as {row, col}
    return shift[1], shift[0]


def dpixels_to_dwcs(dx, dy, pixscale_x, pixscale_y, ref_dec):
    """Converts from pixel shifts to shifts in ra and declination"""
    d_dec = dy * pixscale_y

    d_ra = dx * pixscale_x * np.cos(np.deg2rad(ref_dec))
    return d_ra, d_dec


def get_pixel_scale(header):
    x_scale = header["CD1_1"]
    y_scale = header["CD2_2"]
    return x_scale, y_scale


def get_image_coords(header):
    """Fetches RA and dec from a given header"""
    return header["CRVAL1"], header["CRVAL2"]


