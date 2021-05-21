import glob
import os
import shutil

def make_align_copy(source_files, destination_dir, verbose=False):
    """ Function for making a clean copy
    """
    if verbose:
        print('Copying', source_files)
        print('Destination', destination_dir)

    filelist = []
    for fn in source_files:
        name = os.path.basename(fn)
        out_path = destination_dir + out_dir_name
        if not os.path.isdir(out_path):
            os.mkdir(out_path)
        new_name = out_path + name
        print("copied to")
        print(new_name)
        shutil.copyfile(fn, new_name)
        filelist.append(new_name)

    return filelist

def clean_align_copy(align_directory):
    pass


def write_wcs_to_fits(input_image, dra, ddec):
    """ Writes the new coordinates into all the headers in the fits image.

    Parameters
    ----------
    input_image : str
        path to the image that should be updated

    dra : float
        change to Ra coordinate in decimal degrees
    ddec : float
        change to Dec coordinate in decimal degrees

    Returns
    -------
    None
    """
    print('Propagating changes to: ', input_image,)
    print('Which has the following hdu objects')
    with fits.open(input_image, mode="update") as hdul:
        print(hdul)
        for hdu in hdul:
            try:
                new_ra = hdu.header["CRVAL1"] + dra
                new_dec = hdu.header["CRVAL2"] + ddec
                hdu.header.set("CRVAL1", new_ra)
                hdu.header.set("CRVAL2", new_dec)
                print("Worked: ", hdu)
            except KeyError:
                # This means this HDU did not have a coordinate system
                print(input_image)
                print(hdu)
                continue

        hdul.flush()

