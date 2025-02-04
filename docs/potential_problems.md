# Potential Problems and suggested solutions 


Most problems will appear as insufficient, poor or simply incorrect corrections and unfortunately
this can have many reasons, each of which has it's own set of solutions.

## Window size
The 2D correlation is very sensitive to noise. This means that we want the iage size used for correlation to be as small as possible 
while still containing the main galaxy structure. 

In the case of large offsets one may need to apply manual offsets in order to bring the images einto approximate alignment before 
running the code. This is described [here](https://runnholm.github.io/CROCOA/manual_shifts/).

## Normalization
Before the crosscorrelation is done the images need to be normalized. This can be done in several different ways and certain datasets are more or less 
suitable for a given method, depending on the pixel-statistics of a given individual image.


## CR Rejection
If there are cosmic rays present in an image this will naturally make the cross correlation difficult. This is not an issue for ACS/SBC since the mama
detector does not show cosmic rays. 

CRs need to be masked separately before crosscorrelation and making sure that they are masked properly is crucial. We recommend using `astroscrappy`.

However, when the drizzling is done in CROCOA it is done with a very small pixel scale which can lead to issues where some final pixels only have 
contributions from masked pixels. 




