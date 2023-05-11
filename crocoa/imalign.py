import numpy as np
import pandas as pd
import glob
import os
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from drizzlepac.astrodrizzle import AstroDrizzle as adriz
from crocoa import filemanagement as fm
from astropy.io import fits
from scipy.signal import correlate2d
from skimage.registration import phase_cross_correlation
from datetime import datetime

def align(
    source_images,
    destination_dir,
    run_drizzle=True,
    drizzle_config=None,
    drizzle_groups=None,
    temp_dir="./temp",
    reference_image=None,
    cleanup=True,
    normalization="white",
    hlet=False,
    manual_shift=None,
    mask_method="zeros",
    verbose=False,
):
    """Function for aligning two image sets
    Parameters
    ----------
    source_images : list of str
        List or array containing the images that are to be aligned
    destination_dir : str
        directory in which to put matched files
    run_drizzle : bool
        whether or not to drizzle the frames before comparing. If false the source images are assumed to be drizzled
    drizzle_config : dict
        dictionary with astrodrizzle configuration settings
    drizzle_groups : list of lists
        list containing the groupings of which frames that should be drizzled together.
        If `None`each frame is drizzled separately.
    temp_dir : str
        Directory where intermediate results are stored
    reference_image : str or None
        sets the reference image for the correlation (full path of the image)
    cleanup : bool, default = True
        whether or not to remove the temp_dir
    normalization : str
        Parameter passed to the image_normalization function.
    mask_method : str
        NaN treatment in drizzled frames, 'zero' (default) or 'skimage'.
    hlet : bool
        whether to store the new coord system as an hlet in the fits. Not implemented
    verbose : bool
        whether to show the output of astrodrizzle

    """
    # Step 1: make drizzle copies
    # Create working dir
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)
    else:
        # here we want to deal with an already existing temp directory which
        # can screw us over
        temp_dir = temp_dir + datetime.now().strftime("%m%d%Y_%H%M%S")
        os.mkdir(temp_dir)
        if verbose:
            print('Warning: Previous temp-directory found. Renamin temp to: ' + temp_dir)



    driz_source_dir = temp_dir + "/raw_images/"
    driz_destination_dir = temp_dir + "/drizzled_images/"
    if not os.path.isdir(driz_destination_dir):
        os.mkdir(driz_destination_dir)
    driz_source_files = fm.make_copy(source_images, driz_source_dir, verbose=verbose)

    if manual_shift is not None:
        assert type(manual_shift) is dict
        manual_wcs_shift(driz_source_files, manual_shift)

    # Step 2: drizzle
    if run_drizzle:
        # TODO: generalize this to be able to take groups of files as well
        for fd_frame in driz_source_files:
            if verbose:
                adriz(
                    input=[fd_frame],
                    output=driz_destination_dir
                    + os.path.basename(fd_frame).split(".")[0],
                    **drizzle_config
                )
            else:
                with suppress_stdout_stderr():
                    adriz(
                        input=[fd_frame],
                        output=driz_destination_dir
                        + os.path.basename(fd_frame).split(".")[0],
                        **drizzle_config
                    )

        # Step 4: Copy files to destination dir
        if not os.path.isdir(destination_dir):
            os.mkdir(destination_dir)

    destination_files = fm.make_copy(
        source_images, destination_dir + "/", verbose=verbose
    )

    # Step 3: Cross Correlate
    if run_drizzle:
        corr_source_files = glob.glob(driz_destination_dir + "*sci*")
    else:
        corr_source_files = source_images
    if reference_image is None:
        reference_image = corr_source_files[0]
        corr_source_files = corr_source_files[1:]
    else:
        if run_drizzle:
            ref = glob.glob(driz_destination_dir + os.path.basename(reference_image))
            if ref:
                reference_image = ref[0]
            else:
                print("Cannot find reference image. Falling back to using first image")
                reference_image = corr_source_files[0]
                corr_source_files = corr_source_files[1:]
        else:
            ref = reference_image

    for frame in corr_source_files:
        destination_file = glob.glob(
            destination_dir + "/" + os.path.basename(frame).split("_")[0] + "*"
        )
        match_images(
            reference_image,
            frame,
            destination_file,
            normalization=normalization,
            mask_method=mask_method,
        )
        if manual_shift is not None:
            print()
            print("manual shift")
            manual_wcs_shift(destination_file, manual_shift)
            print()

    # Add the manual shift to the reference file
    if manual_shift is not None:
        reference_destination = destination_file = glob.glob(
            destination_dir
            + "/"
            + os.path.basename(reference_image).split("_")[0]
            + "*"
        )
        manual_wcs_shift(reference_destination, manual_shift)

    # Step 5: Cleanup
    if cleanup:
        fm.cleanup(temp_dir)


def match_images(
    reference_fits,
    image_fits,
    output_fits,
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
    output_fits : list
        path to original flc files where wcs should be updated
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

    # 4. Write new coordinates to output
    for outfile in output_fits:
        if verbose:
            print("writing to: ", outfile)
        fm.write_wcs_to_fits(outfile, dra, ddec)


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


@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(os.devnull, "w") as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)
