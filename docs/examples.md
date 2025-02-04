# Examples
Here are some examples of primary use cases for the code


The examples here assumes that the flt frames have been pipeline processed and 

## Aligning frames within a single filter

Add the central Ra and Dec to the drizzle configuration 

``` py
from crocoa.drizzle_config import config

config['final_ra'] = 1.12637
config['final_dec'] = -10.19156
```

Set up the working directory and create an `ImageSet`

!!! warning
    This will change the world coordinate systems of the flt files that are given to the ImageSet. 
    Making copies is always advisable

``` py
from crocoa.filemanagement import ImageSet
from pathlib import Path

source_dir = Path('path/to/flts')
flt_frames = list(source_dir.glob('*.flt'))
img_set = ImageSet(flt_frames, filtername="F150LP", drizzle_config=config)
```

Use the align filter function to align the images in the `ImageSet`
``` py
from crocoa.imalign import align_single_filter

align_single_filter(img_set)

```


## Aligning frames between several filters
This is done in a very similar manner to aligning a single filter, however the difference is that 
we pass several ImageSets as a list rather than a single image set

!!! note
    This assumes that each individual filter has already been aligned internally according to the process
    outlined above.

``` py
from crocoa.filemanagement import ImageSet
from crocoa.imalign import align_multiple_filters
from crocoa.drizzle_config import config

from pathlib import Path

config['final_ra'] = 1.12637
config['final_dec'] = -10.19156

# Filter 1
source_dir = Path('path/to/filter1/flts')
flt_frames = list(source_dir.glob('*.flt'))
img_set1 = ImageSet(flt_frames, filtername="filter1", drizzle_config=config)

# Filter 2
source_dir = Path('path/to/filter2/flts')
flt_frames = list(source_dir.glob('*.flt'))
img_set2 = ImageSet(flt_frames, filtername="filter2", drizzle_config=config)

img_sets = [img_set1, img_set2]

align_multiple_filters(img_sets)
```

## Tweaking configuration
If you want to pass more detailed configurations to the internal steps this is done in 2 ways

### Changing Drizzle configuration settings
The default drizzle configuration settings used in the code are located in `crocoa.drizzle_config`
and are passed as a dictionary

If you wish to examine or alter these you can simply import these settings, alter them and pass the
new dictionary to the ImageSet constructor when you set it up

``` py
from crocoa.filemanagement import ImageSet
from crocoa.drizzle_config import config

# Examine the values in the current config
for k in config.keys():
    print(k, ': ', config[k])

# Change one 
config[final_scale] = 0.02

# Then pass it to your image set 
```

### Changing image correlation settings
The settings available for the image correlation are the keywords specified for the `match_images` function

- `verbose=True`
- `normalization="white"`
- `clip=None`
- `mask_method="zero"`

The verbosity sets the detail of the printout for each matching. When not debugging it can be useful to set this to false to prevent 
an enormous amount of terminal output from astro_drizzle

The normalization sets what method is used to normalize the images before running the cross correlation. 

The clip argument sets whether some range of values should be clipped out of the image. 

The mask method sets what way nan values in the resulting drizzled images should be handled. By default they are just set to zero


To change these arguments one passes a dictionary to the align function 

``` py
from crocoa.imalign import align_single_filter

# Lets use a minmax normalization in the matching 
crosscorr_config = {'normalization': 'minmax'}

align_single_filter(img_set, matching_config=crosscorr_config)

```
