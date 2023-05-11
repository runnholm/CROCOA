import glob
import os
import shutil
from astropy.io import fits


def make_copy(source_files, destination_dir, verbose=False):
    """Function for making a clean copy"""
    if verbose:
        print("Copying", source_files)
        print("Destination", destination_dir)

    filelist = []
    for fn in source_files:
        name = os.path.basename(fn)
        if not os.path.isdir(destination_dir):
            os.mkdir(destination_dir)
        new_name = destination_dir + name
        print("copied to")
        print(new_name)
        shutil.copyfile(fn, new_name)
        filelist.append(new_name)

    return filelist


def cleanup(directory):
    shutil.rmtree(directory)


def write_wcs_to_fits(input_image, dra, ddec):
    """Writes the new coordinates into all the headers in the fits image.

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
    print(
        "Propagating changes to: ",
        input_image,
    )
    print("Which has the following hdu objects")
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


class Image:
    def __init__(self, filename, verbose=False) -> None:
        self.original = filename
        self.verbose = verbose
    
    def make_working_copy(self, destination_dir):
        filelist = make_copy([self.original], destination_dir, verbose=self.verbose)
        self.working_copy = filelist[0]
    
    def backpropagate_wcs(self, dra, ddec):
        write_wcs_to_fits(self.original, dra, ddec)
    

class ImageSet:
    def __init__(self, file_list, drizzle_config) -> None:
        self.images = [Image(file_name) for file_name in file_list]
        self.drizz_result = None
        self.drizz_config = drizzle_config


    def drizzle(self):
        pass


def align_multiple_filters(image_sets, reference_set_index=0):
    for image_set in image_sets:
        image_set.drizzle()
    
    



    