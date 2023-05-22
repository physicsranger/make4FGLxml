# make4FGLxml
An updated version of the make4FGLxml script, including a GUI, used to produce source region models for analysis of Fermi LAT data.  If you have previous experience with the script (pre-GitHub), consider this version 1.10.0 (or v01r10 using my old syntax).

The main changes are to improve the readability of the code and make it more elegant (though elegence is a matter of opinion).  If you are already familiar with the code, the [previous instructions](https://fermi.gsfc.nasa.gov/ssc/data/analysis/user/readme_make4FGLxml.txt) are still mostly valid, with the major changes being in some argument names and that the custom class is called SourceList instead of srcList.  Additionally, previously the sources were grouped by distance from the ROI center with comments inserted to denote the "distance block", but only when using the FITS catalog.  Currently, the script no longer does this, but it does sort the sources by distance from the ROI center, even if you're using the XML version of the catalog, which it didn't do previously.

The script now requires the ```pandas``` module, which can easily be installed with pip or your preferred package manager.

---

## Installation
Not much to say here yet.  Clone the repo and make sure it is in your PYTHONPATH and PATH if you want to use the script as a command line executable. Note that, in addition to the main make4FGLxml.py script, there are two submodules, build_model and make4FGLxml_Components, which need to be in your python path so simply copying the main script to a directory in your PYTHONPATH won't work.

I will eventually include instructions for installing the GUI as a stand alone executable (using pyInstaller) and I might arrange things to be installable via pip.

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

```python
#using an event file to get the ROI information and the FITS version of the catalog
#and not specifying either optional parameter
source_list=SourceList(catalog_file='/some/path/to/gll_psc_v31.fit','my_event_file.fits')

#or...

#using a list with ROI information adn the XML version of the catalog
source_list=SourceList(catalog_file='/some/path/to/gll_psc_v28.xml',[123.4,-32.2,15],
                       output_name='my_LAT_model.xml',
                       write_directory='/some/path/to/analysis/directory')
```
Once this object is created, the only thing left do is call the ```make_model``` method.  This method has several arguments, but all of them are optional, though some manual editing of the model may be needed if the default parameters are used.  The ```make_model``` arguments are:
 * _galactic\_file_ (str, optional): path to the Galactic diffuse model template, if not specified the code assumes gll\_iem\_v07.fits and you will have to add the full path information
 * _galactic\_name_ (str, optional): name of the Galactic diffuse component in the model
 * _isotropic\_file_ (str, optional): path to the isotropic component template, it not specified the code assumed iso\_P8R3\_SOURCE\_V3\_v1.txt and you will have to add the full path information
 * _isotropic\_name_ (str, optional): name of the isotropic diffuse component in the model
 * _norms\_free\_only_ (bool, optional): flag to allow only normalization parameters to be free for sources which satisfy other requirements to be free, default is False
 * _extended\_directory_ (str, optional): path to directory with extended source templates, if not specified the code will assume the default location within the _fermitools_
 * _free\_radius_ (float, optional): radial distance from ROI center within which to free parameters of sources which meet significance requirements, if not specified will default to ROI radius
 * _max\_free\_radius_ (float, optional): maximum radial distance From ROI center for any free sources to have free parameters, only applicable if freeing variable sources
 * _extra\_radius_ (float, optional): radial extent beyond ROI radius to include sources, with all parameters fixed, to account for the large LAT PSF at low energies, defaults to 10 degrees
 * _sigma\_to\_free_ (float, optional): minimum average significance (test statistic), if using the FITS (XML) catalog, a source must have to consider setting parameters free
 * _variable\_free_ (bool, optional): flag to free the normalization parameters of sources found to be significantly variable, even if the source does not meet significance requirements or is > _free\_radius_ degrees but $\leq$ _max\_free\_radius_ away from the ROI center, default is True
 * _force\_point\_sources_ (bool, optional): flag to include extended sources as point sources in the model, default is False
 * _extended\_catalog\_names_ (bool, optional): flag to use the "4FGL JXXXX.X+XXXXe" names for extended sources, default is False, note that this only applies when using the FITS version of the catalog, if you use teh XML version then the behavior is always as if this flag was set to True
 * _make\_region_ (bool, optional): flag to also generate a ds9 region file, default is True
 * _region\_file_ (str, optional): name of output ds9 region file, will be written in _write\_directory_, if not specified it will prepend 'ROI_' to _output\_name_ and change the extension to '.reg'
 * _galactic\_index\_free_ (bool, optional): flag to modify the spectrum of the Galactic diffuse emission using a power-law model with free Index parameter, default is True
 * _use\_old\_names_ (bool, optional): flag to use source names following the convention of make1FGLxml.py and make2FGLxml.py, i.e., "\_4FGLJXXXX.X+XXXX", default is False

In general, most users will be able to use the defaults for all of the parameters except _free\_radius_, _max\_free\_radius_, and _sigma\_to\_free_.  The example below shows the parameter values I typically use when preparing for a binned likelihood analysis using an ROI with 15 degree radius and the FITS version of the catalog:

```python
source_list.make_model(free_radius=6,max_free_radius=8,sigma_to_free=12)
```

As noted, the previous example assumes that the ```source_list``` object was created referencing the FITS version of the 4FGL catalog, so the values given mean that sources only have free spectral parameters if they were found in the catalog with $\geq12\sigma$ average significance and are within 6 degrees of the ROI center or if they were found to be significantly variable in the catalog and are within 8 degrees of the ROI center.  If the ```source_list``` object had been created referencing the XML version of the catalog, the _sigma\_to\_free_ parameter applies to the source test statistic, so a value more like 100 or 200 might be more suitable than 12, depending on the use case.

### Using the Command Line Interface
The script can be called as an executable from the command line.  You'll need to add the repo to your PATH to be able to call it.  You can then invoke the script either on it's own:

```make4FGLxml.py --options```  

or using python:

```python make4FGLxml.py --options```  

If the first option doesn't work, you should try editing the shebang at the top of the file to point to the correct python on your system.

The command line options are almost entirely the same, with the same names, as those described in the section for how to use the script in an interactive session or your own script.  The only difference is in how the ROI information is entered. To get the ROI information from a _Fermi_ LAT event file (in this case we'll assume it is called my\_LAT\_events.fits and is in the current directory), use a command similar to:

```make4FGLxml.py /some/path/to/gll_psc_v28.xml --event_file my_LAT_events.fits --options```

To instead provide the ROI information directly (here we will center at (RA,DEC) = (120.3 deg., 30.5 deg.) with a 15 degree radius), use a command similar to:

```make4FGLxml.py /some/path/to/gll_psc_v31.fits --RA 120.3 --DEC 30.5 --radius 15 --options```

When specifying optional arguments (represented by ```--options``` in the example commands above)  you can use the same names as in the interactive example, preceded by two dashes.  Some options do have shortened identifiers preceded by only one dash.  To see all the options and help text, simply use the command:

```make4FGLxml.py --help```

### Using the GUI
The GUI is still a work in progress and has not been fully tested; however, early testing does suggest it is mostly working, but use with caution.  Currently, you can launch the GUI from an interactive python session:

```python
from make4FGLxml_GUI import main
main()
```

I am currently working on getting things to actually work from the command line, and eventually for the ability to make a standalone executable with pyInstaller.
