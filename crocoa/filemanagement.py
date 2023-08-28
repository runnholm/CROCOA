import glob
import os
import shutil
from astropy.io import fits
from pathlib import Path
from drizzlepac.astrodrizzle import AstroDrizzle as adriz
from crocoa.utilities import suppress_stdout_stderr
import glob


def make_copy(source_files, destination_dir, target_name=None, verbose=False):
    """Function for making a clean copy"""
    if verbose:
        print("Copying", source_files)
        print("Destination", destination_dir)

    filelist = []
    for fn in source_files:
        if target_name is None:
            name = os.path.basename(fn)
        else:
            name = target_name
        if not os.path.isdir(destination_dir):
            os.mkdir(destination_dir)
        if isinstance(destination_dir, Path):
            new_name = destination_dir / name
        else:
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
    def __init__(self, filename, verbose=False, file_suffix=None) -> None:
        basename = os.path.basename(filename)
        if file_suffix is None:
            self.name = basename
        else:
            # put suffix on the right side of flt
            if 'flt' in basename:
                self.name = basename.split("flt")[0] + file_suffix + '_flt.' + basename.split(".")[1]
            if 'flc' in basename:
                self.name = basename.split("flc")[0] + file_suffix + '_flc.' + basename.split(".")[1]
            else:
                self.name = basename.split(".")[0] + file_suffix + basename.split(".")[1]
        self.original = filename
        self.verbose = verbose
        self.target_copy = None

    def make_working_copy(self, destination_dir):
        filelist = make_copy(
            [self.original],
            destination_dir,
            target_name=self.name,
            verbose=self.verbose,
        )
        self.working_copy = filelist[0]

    def make_target_copy(self, destination_dir):
        filelist = make_copy(
            [self.original],
            destination_dir,
            target_name=self.name,
            verbose=self.verbose,
        )
        self.target_copy = filelist[0]

    def backpropagate_wcs(self, dra, ddec):
        if self.target_copy is not None:
            write_wcs_to_fits(self.target_copy, dra, ddec)
        else:
            write_wcs_to_fits(self.original, dra, ddec)

    def shift_working_copy_wcs(self, dra, ddec):
        """Utility function for modifying the working copy wcs"""
        write_wcs_to_fits(self.working_copy, dra, ddec)


class ImageSet:
    def __init__(
        self,
        file_list,
        filtername,
        drizzle_config,
        working_dir="./temp",
        destination_dir=None,
        manual_shifts=None,
        verbose=False,
        file_suffix=None,
    ) -> None:
        """
        Parameters
        ----------
        file_list : list
            list of string file names to be included in the image set
        filtername : str
            Name of the current filter - just used to create subdirectories
        drizzle_config : dict
            configuration dictionary for astrodrizzle
        working_dir : str, optional
            sets the temp directory parent name
        destination_dir : str, optional
            This sets the directory where the resulting files with the backpropagated wcs' are put. If None this defaults
            to './filtername_aligned'
        manual_shifts : dict or tuple
            Either a dictionary mapping dra ddec shifts to individual images in the set or a tuple with a global dra and
            ddec that should be applied to all the frames collectively
        verbose : bool
            whether to output drizzle resutls
        file_suffix : str, optional
            suffix to append to the file name(s)
        """
        self.images = [
            Image(file_name, verbose=verbose, file_suffix=file_suffix)
            for file_name in file_list
        ]
        self.verbose = verbose
        self.drz_config = drizzle_config
        self.filtername = filtername
        self.drz_source_dir = Path(working_dir) / "drz_source" / filtername
        self.drz_target_dir = Path(working_dir) / "drz_target" / filtername
        self.manual_shifts = manual_shifts

        # Set up destination files
        self.destination_dir = destination_dir
        if self.destination_dir is not None:
            self.destination_dir = Path(destination_dir)
            self.destination_dir.mkdir(parents=True, exist_ok=True)
        else:
            self.destination_dir = Path("./{}_aligned".format(filter))
            self.destination_dir.mkdir(parents=True, exist_ok=True)

        # set up drizzle files
        self.drz_source_dir.mkdir(parents=True, exist_ok=True)
        self.drz_target_dir.mkdir(parents=True, exist_ok=True)

    def make_all_copies(self):
        """Create all relevant file copies"""
        self.make_target_copies()
        self.make_driz_source()

    def make_target_copies(self):
        """Make target file copies"""
        for img in self.images:
            img.make_target_copy(self.destination_dir)

    def make_driz_source(self):
        """Make source file copies since drizzle changes the source files it's run on"""
        self.working_source = []
        for image in self.images:
            image.make_working_copy(self.drz_source_dir)
            self.working_source.append(str(image.working_copy))

    def drizzle(self, individual=False):
        """Run the drizzling process"""
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
                                self.drz_target_dir
                                / os.path.basename(image).split(".")[0]
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

    def apply_manual_shifts(self):
        """Take a manual shift and apply it to either all images collectively or to each individual frame"""
        # First parse to understand if we have a collective change or not
        if self.manual_shifts is None:
            if self.verbose:
                print("No manual shifts specified")
            pass

        if isinstance(self.manual_shifts, dict):
            for image in self.images:
                if image in self.manual_shifts.keys():
                    # Add the shift to the working copy
                    image.shift_working_copy_wcs(
                        self.manual_shifts[image]["dra"],
                        self.manual_shifts[image]["ddec"],
                    )
                    # also add it to the target
                    image.backpropagate_wcs(
                        self.manual_shifts[image]["dra"],
                        self.manual_shifts[image]["ddec"],
                    )
        elif isinstance(self.manual_shifts, tuple):
            for image in self.images:
                image.shift_working_copy_wcs(
                    self.manual_shifts[0], self.manual_shifts[1]
                )
                image.backpropagate_wcs(self.manual_shifts[0], self.manual_shifts[1])

    def backpropagate_wcs_shift(self, dra, ddec):
        """propagate the found wcs offset to all images"""
        for image in self.images:
            image.backpropagate_wcs(dra, ddec)

    def clean_temp_directories(self):
        """Remove drizzling directories"""
        shutil.rmtree(self.drz_source_dir)
        shutil.rmtree(self.drz_target_dir)

    def __repr__(self) -> str:
        msg = "< {} ImageSet instance containing: ".format(self.filtername)
        for img in self.images:
            msg += img + ", "
        msg += ">"
        return msg
