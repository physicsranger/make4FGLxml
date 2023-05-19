# make4FGLxml
An updated version of the make4FGLxml script, including a GUI, used to produce source region models for analysis of Fermi LAT data.

The main changes are to improve the readability of the code and make it more elegant (though elegence is a matter of opinion).  If you 
are already familiar with the code, the [previous instructions](https://fermi.gsfc.nasa.gov/ssc/data/analysis/user/readme_make4FGLxml.txt) are still mostly valid, with the major
changes being in some argument names and that the custom class is called SourceList instead of srcList.

---

## Installation
Not much to say here yet.  Clone the repo and make sure it is in your PYTHONPATH and PATH if you want to use the script as a command line executable.
Note that, in addition to the main make4FGLxml.py script, there are two submodules, build_model and make4FGLxml_Components, which need to be
in your python path so simply copying the main script to a directory in your PYTHONPATH won't work.

I will eventually include instructions for installing the GUI as a stand alone executable (using pyInstaller) and I might arrange things to be
installable via pip.

---

## Using the Script...
### In an Interactive Python Session or Personal Script
To use the script, simply import the custom ```SourceList``` class in an interactive session or your own script as shown below.

```python
from make4FGLxml import SourceList
```

When creating a ```SourceList``` object, there are two required arguments and two optional arguments:
 * _catalog\_file_ (str): path to a 4FGL XML or FITS catalog file
 * _ROI_ (str or list): either the path to a _Fermi_ LAT event file for the region you plan to analyze, or a list with the RA, DEC, and radius of the region the model should be built for
 * _output\_name_ (str, optional): the name of the resulting spatial-spectral XML model
 * _write\_directory_ (str, optional): the path to the desired output directory for _output\_name_.

See the two example usages below:
