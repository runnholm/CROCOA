# Adding manual shifts to the images

In many cases the wcs of a filter or even an individual frame may be significantly offset from the
others in a set. In this case one may be forced to use a relatively large region in the drizzling
configuration. To avoid this the code allows for adding shifts to the RA and dec manually both
to entire filters and to individual frames.



## Shifting individual frames
It is not uncommon that a given exposure may be significantly shifted from the others in a set 
especially if exposures were taken in different visits. In this case the `manual_shifts` parameter
needs to be a dictionary specifying the full filename (not full path) of the frame that needs shifting.
for each such frame the shift is then given by a nested dict with the `'dra'` and `'ddec'` keys. 

The offsets are given in decimal degrees

``` py
manual_shifts = {
    'full_filename_flt.fits':{
        'dra': 0.02,
        'ddec': 0.01,
    }, 
    'full_filename2_flt.fits':{
        'dra': -0.02,
        'ddec': 0.01,
    }
}

```
One does not have to specify shifts for all of the frames in a given image set. This is then given to the 
`ImageSet`during construction


## Adding a shift to allframes in a filter
This is only relevant after the frames in the filter have been aligned internally but instead of giving the
shifts for individual frame names, they are simply specified as tuple of `(dra, ddec)` given to the relevant ImageSet


``` py

from crocoa.filemanagement import ImageSet
from crocoa.imalign import align_multiple_filters
from crocoa.drizzle_config import config

from pathlib import Path

config['final_ra'] = 1.12637
config['final_dec'] = -10.19156
manual_shifts = (0.002, 0.2)


# Filter 1
source_dir = Path('path/to/filter1/flts')
flt_frames = list(source_dir.glob('*.flt'))
img_set1 = ImageSet(flt_frames, drizzle_config=config, manual_shifts=manual_shifts)


... 
# Remaining filters etc
``` 



