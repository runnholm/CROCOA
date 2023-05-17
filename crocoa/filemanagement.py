import glob
import os
import shutil
from astropy.io import fits
from pathlib import Path
from drizzlepac.astrodrizzle import AstroDrizzle as adriz
from crocoa.utilities import suppress_stdout_stderr
import glob


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
    def __init__(
        self, filtername, file_list, drizzle_config, working_dir="./temp", verbose=False
    ) -> None:
        """
        filtername : str
            Name of the current filter - just used to create subdirectories
        file_list : list
            list of string file names to be included in the image set
        drizzle_config : dict
            configuration dictionary for astrodrizzle
        working_dir : str, optional
            sets the temp directory parent name
        verbose : bool
            whether to output drizzle resutls
        """
        self.images = [Image(file_name, verbose=verbose) for file_name in file_list]
        self.verbose = verbose
        self.drz_config = drizzle_config
        self.filtername = filtername
        self.drz_source_dir = Path(working_dir) / "drz_source" / filtername
        self.drz_target_dir = Path(working_dir) / "drz_target" / filtername

        self.drz_source_dir.mkdir(parents=True)
        self.drz_target_dir.mkdir(parents=True)

    def make_driz_source(self):
        self.working_source = []
        for image in self.images:
            image.make_working_copy(self.drz_source_dir)
            self.working_source.append(image.working_copy)

    def drizzle(self, individual=False):
        self.drizzled_files = []
        if individual:
            for image in self.working_source:
                if self.verbose:
                    adriz(
                        input=[image],
                        output=str(
                            self.drz_target_dir / os.path.basename(image).split(".")[0]
                        ),
                        **self.drz_config
                    )
                else:
                    with suppress_stdout_stderr():
                        adriz(
                            input=[image],
                            output=str(
                                self.drz_target_dir / os.path.basename(image).split(".")[0]
                            ),
                            **self.drz_config
                        )
                self.drizzled_files.append(
                    self.drz_target_dir.glob(
                        os.path.basename(image).split(".")[0] + "*drz*.fits"
                    )
                )
        else:
            if self.verbose:
                adriz(
                    input=self.working_source,
                    output=str(self.drz_target_dir / self.filtername),
                    **self.drz_config
                )
            else:
                with suppress_stdout_stderr():
                    adriz(
                        input=self.working_source,
                        output=str(self.drz_target_dir / self.filtername),
                        **self.drz_config
                    )
            # get the resulting drizzled file(s)
            self.drizzled_files = self.drz_target_dir.glob("*drz*.fits")

    def backpropagate_wcs_shift(self, dra, ddec):
        for image in self.images:
            image.backpropagate_wcs(dra, ddec)

    def clean_temp_directories(self):
        shutil.rmtree(self.drz_source_dir)
        shutil.rmtree(self.drz_target_dir)



