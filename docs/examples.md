# Examples
Here are some examples of primary use cases for the code


The examples here assumes that the flt frames have been pipeline processed and 

## Aligning frames within a single filter

1. Set up the working directory and create an `ImageSet`
``` py
from crocoa.filemanagement import ImageSet
from pathlib import Path

source_dir = Path('path/to/flts')
flt_frames = list(source_dir.glob('*.flt'))
img_set = ImageSet(flt_frames)
```
2. Use the align filter function to align the images in the `ImageSet`
``` py
from crocoa.imalign import align_single_filter

align_single_filter(img_set)

```




## Aligning frames between several filters

