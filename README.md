# CROCOA
<b>CRoss COrrelation Alignment tool. </b>

Software package for using 2D crosscorrelation to align astronomical images.



## Requirements
The primary requirements are astropy and scipy in order to manage the fits files and 2d correlation function. Additionally drizzlepac is also used to drizzle the frames in order to account for differing pixel grids and rotation

The package is tested on python 3.7


## Installation

For having an editable installation activate your preferred virtualenv and run the following in the CROCOA main directory

``` 
$ python setup.py develop
```

## Usage
Example usage:

```python
from crocoa.imalign import align
from crocoa.drizzle_config import config

# Review the drizzle config settings
for k in config.keywords():
    print(k, ': ', config[k])

# Set the Ra and dec in the drizzle config to accurate values
config['final_ra'] = 1.12637
config['final_dec'] = -10.19156

# Prepare your directories
source_images = glob.glob("SDSSJ0004-1011/F150LP/*")

destination_dir = "SDSSJ0004-1011/F150LP_aligned/"
os.mkdir(destination_dir)

# Run alignment
align(source_images, destination_dir, config)
```
## Contributions

Contributions are welcome in the form of pull requests and issue submissions.