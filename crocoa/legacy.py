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

from crocoa.imalign import *

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
            print('Warning: Previous temp-directory found. Renaming temp to: ' + temp_dir)



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
