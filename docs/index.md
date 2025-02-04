# CROCOA
<b>CRoss COrrelation Alignment tool. </b>

Software package for using 2D crosscorrelation to align astronomical images.



## Requirements
The primary requirements are astropy and scipy in order to manage the fits files and 2d correlation function. Additionally drizzlepac is also used to drizzle the frames in order to account for differing pixel grids and rotation

The package is tested on python 3.7 and 3.10


## Installation
begin by cloning the code into a folder of your choice 
For having an editable installation activate your preferred virtualenv and run the following in the CROCOA main directory

``` 
$ python setup.py develop
```

## Usage
 
See the [examples page](https://runnholm.github.io/CROCOA/examples/) for concrete examples of how to use the code.




### Tips and tricks

#### Image size
The 2D correlation is very sensitive to the relative size of the target to the frame i.e. a small target in a large frame can give poor results. Therefore it can be useful to minimize the size of the frame considered while still makeing sure that the full target is contained.

In cases where the alignment is very poor this can require manual correction of the coordinate system beforehand.

#### Image normalization
How the image is normalized before being run through the 2D correlation can also impact the quality of the alignment. test setting the normalization keyword to different settings.

#### Always check manually
What the title says. The system is most definitely not infallible.
