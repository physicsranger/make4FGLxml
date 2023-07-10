# LATSourceModel
A module with classes and utility functions to build spatial-spectral XML model files for analysis of data from the _Fermi Gamma-ray Space Telescope_ Large Area Telescope.  The models are constructed from a FITS or XML version of the 4FGL catalog (data release versions 1, 2, 3, & 4).

Development happens on the [make4FGLxml GitHub](https://github.com/physicsranger/make4FGLxml) page.  The module is based on the original make4FGLxml.py script (versions <=1.09.0, and versions for previous LAT catalogs) released as a user contribution through the [Fermi Science Support Center](https://fermi.gsfc.nasa.gov/ssc/data/analysis/user/) user contributions page.

This module updates and streamlines the code while adding new functionality.  Versions 1.10.0 of the make4FGLxml.py script are now only a way to quickly use the new module to make a source model via the command line.

## Authors
 - [@Tyrel Johnson](https://github.com/physicsranger)
 
## Installation

Install with pip

```bash
pip install LATSourceModel
```
## Requirements
* numpy
* pandas
* astropy 