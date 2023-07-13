import astropy.io.fits as pyfits

import numpy as np

from xml.dom import minidom

from pathlib import Path

import warnings

#a list of the available spectral models
#note, this includes all listed on the FSSC
#but the more common models are 'up front'
spectral_models=['PowerLaw','LogParabola','PLSuperExpCutoff4',
                 'PowerLaw2','BrokenPowerLaw','PLSuperExpCutoff2',
                 'FileFunction','BrokenPowerLaw2','SmoothBrokenPowerLaw',
                 'PLSuperExpCutoff','PLSuperExpCutoff3','ExpCutoff',
                 'BPLExpCutoff','Gaussian','ConstantValue',
                 'BandFunction','DMFitFunction']

#a list of available spatial models
spatial_models=['SkyDir','RadialDisk','RadialGaussian','SpatialMap',
                'ConstantValue','MapCubeFunction']

#copied from Damien's macro
#Create an XML document parameter description element
def parameter_element(free,name,maximum,minimum,scale,value):
    '''
    function to create a XML parameter element for a spatialModel
    our spectrum element of a LAT XML model
    
    this function represents slight modifications of
    a function written by Damien Parent while at NRL

    all parameters must eventually be converted to strings, so
    there is flexibility in how they are input

    Parameters
    ----------
    free - int or str
        flag to indicate if parameter should be left free to vary, 1 = True
        and 0 = False
    name - str
        name of parameter
    maximum - float or str
        maximum value parameter can take
    minimum - float or str
        minimum value parameter can take
    scale - float or str
        scale of parameter, parameter is evaluated as
        scale * value
    value - float or str
        value of parameter, after accounting for scale parameter

    Returns
    -------
    XML document element
        XML document representation of the parameter, to be added
        to a spatialModel or spectrum document element
    '''
    
    #get 'implementation'
    impl=minidom.getDOMImplementation()

    #make a 'document
    xmldoc_out=impl.createDocument(None,None,None)

    #create the element
    parameter=xmldoc_out.createElement('parameter')

    #set the element parameters
    parameter.setAttribute('free',str(free))
    parameter.setAttribute('name',str(name))
    parameter.setAttribute('max',str(maximum))
    parameter.setAttribute('min',str(minimum))
    parameter.setAttribute('scale',str(scale))
    parameter.setAttribute('value',str(value))
    
    return parameter

def get_ROI_from_event_file(event_file_name):
    '''
    function to extract the region of interest information
    from a Fermi LAT event FITS file

    Parameters
    ----------
    event_file_name : str or path-like object
        name of the event file from which ROI information is desired
        if the full path is not specified it must be in the current
        directory

    Returns
    -------
    list
        list of floats corresponding to right ascension, declination,
        and radius (all in degrees) from the event file header, or a
        list of None-types if the information was not found in the header
    '''
    
    #open the file and get the DSKEY info from the header
    with pyfits.open(event_file_name) as event_file:
        header=event_file[1].header
        num=header['NDSKEYS']
                
        dstypes={header[f'DSTYP{idx+1}']: header[f'DSVAL{idx+1}']\
            for idx in range(num)}
    
    #now attempt to get the position information and return
    try:
        RA,DEC,radius=dstypes.pop('POS(RA,DEC)').lower().strip('circle()').split(',')
        return float(RA),float(DEC),float(radius)
    
    except KeyError:
        print(f"No position keyword found in header of file {event_file_name}\
 (assuming position is RA and DEC).")
        return None,None,None


#function to calculate the angular separation between a reference location
#on the sky and either single location on the sky or an array of sky locations
#takes arguments in degrees, returns values in degrees
def angular_separation(ref_RA,ref_DEC,RA_values,DEC_values):
    '''
    function to calculate the angular separation between a reference
    position on the sky and either a different, single position
    on the sky or an array of sky locations

    Parameters
    ----------
    ref_RA - float
        right ascension, in degrees, of the reference sky position
    ref_DEC - float
        declination, in degrees, of the reference sky position
    RA_values - array or float
        right ascension value(s) of comparison sky position(s)
    DEC_values - array or float
        declination value(s) of copmarison sky position(s)

    Returns
    -------
    array
        array of angular separations, in degrees, if
        RA_values and DEC_values are single values then
        a length 1 array is returned
    '''
    
    #for ease, create a degree-to-radians conversion variable
    d2r=np.pi/180
    
    #round the number from the calculation to avoid issues with the
    #arccos function
    return np.arccos((np.cos(ref_DEC*d2r)\
        *np.cos(DEC_values*d2r)*np.cos((ref_RA-RA_values)*d2r)\
        +np.sin(ref_DEC*d2r)*np.sin(DEC_values*d2r)).round(6))/d2r
    
#function to create a ds9 region file
def build_region(region_file_name,sources):
    '''
    function to create a ds9-style .reg file from a dictionary
    of sources

    Parameters
    ----------
    region_file_name - str or Path-like
        name of output file, if not the full path then it will
        be created in the current directory
    sources - dict
        nested dictionary of sources from a SourceList class object
    '''
    
    
    #build a shape dictionary and a color list
    shape={'PLSuperExpCutoff':'cross','LogParabola':'diamond',
        'PowerLaw':'circle'}
    color=('magenta','green')

    region_lines='# Region file format: DS9\n'
    region_lines+='# Created by make4FGLxml\n'
    region_lines+='global font="roman 10 normal" move =0\n'

    #cycle through the sources dictionary
    for source_name in sources.keys():
        source=sources.get(source_name)
        entry=f"J2000;point({source['RA']:.3f},{source['DEC']:.3f}) # point = "
            
        #check if source is extended or not
        if source['Extended']:
            entry+=f"box 18 color = {color[source['free']]} "
            
        else:
            entry+=f"{shape[source['spectrum']['model'].strip('2345')]} 15 "
            entry+=f"color = {color[source['free']]} "
            
        #add the name and write to file
        entry+=f"text={{{source_name}}}\n"
        region_lines+=entry

    Path(region_file_name).write_text(region_lines)

def add_to_reg_file(region_file,source_name,source_spatial,spectral_model,extended=False,new_region=None):
    '''
    function to add a new source to an existing region file and either update the file
    or create a new file

    Paramters
    ---------
    region_file - str or Path-like
        name of or path to ds9-style .reg file to be updated
    source_name - str
        name of source to be added
    source_spatial - dict
        dictionary with spatial information, namely RA and DEC which should be the
        right ascension and declination (J2000) of the source to be added
    spectral_model - str
        name of the spectral model used for the source to be added
    extended - bool
        flag to indicate if the source to be added is an extended source (True) or
        a point source (False)
    new_region - str, Path-like, or None-type
        name of or path to the new ds9-style .reg file to be created, if None-type then
        the existing region_file will be overwritten
    '''

    #check that the spatial info has RA and DEC
    if not {'RA','DEC'}.issubset(source_spatial.keys()):
        warnings.warn('Input source spatial information does not have RA and/orDEC.\
 Not adding source to region file.')

    #check if the file passed in actually exists
    elif not Path(region_file).exists():
        warnings.warn(f'Input {str(region_file)} does not exists.')

    else:
        new_line=f"J2000;point({source_spatial['RA']:.3f},{source_spatial['DEC']:.3f}) # point = "

        #follow the point style and size conventions, allowing for spectral models not
        #used in 4FGL, and set the color to cyan to distinguish non-4FGL sources
        if extended:
            new_line+=f"box 18 color = cyan "

        else:
            shape='cross' if spectral_model.strip('2345')=='PLSuperExponentialCutoff' else \
                  'diamond' if spectral_model=='LogParabola' else \
                  'circle' if spectral_model.strip('2')=='PowerLaw' else \
                  'x'

            new_line+=f"{shape} 15 color = cyan "

        new_line+=f"text = {{{source_name}}}\n"
        
        #cet the current .reg file lines
        region_lines=Path(region_file).read_text(encoding='UTF-8')

        #add the new source at the end (each line ends with a newline character)
        region_lines+=new_line

        #check if new_region is None-type or not
        new_region=new_region if new_region is not None else region_file

        #write out the .reg info with the new source added
        Path(new_region).write_text(region_lines)

def valid_spatial_input(spatial_info):
    '''
    function to check if the input dict has valid parameter
    names and all necessary parameters for a Spatial class object

    Parameters
    ----------
    spatial_info - dict
        dictionary with keys and values to create a
        Spatial class object (see build_mode.model_components
        module for more inforamtion)

    Returns
    -------
    bool
        True if all parameter names are valid and all necessary
        parameters are specified, otherwise False
    '''
    
    #check if the model requested is valid
    model=spatial_info.get('spatial_model')
    valid_model=model in spatial_models

    #check if necessary parameters are specified
    #NOTE: this does NOT check the values supplied

    #first, do a check on thoe models with RA and DEC
    if model in ['SkyDir','RadialDisk','RadialGaussian']:
        valid_pars=spatial_info.get('RA') is not None and \
                    spatial_info.get('DEC') is not None

        #check on size parameters for Radial models
        if model!='SkyDir':
            #if the previous logic was True
            #check that one of Radius or Sigma is specified
            valid_pars=valid_pars and \
                        (spatial_info.get('Radius') is not None or \
                         spatial_info.get('Sigma') is not None)

    #next, check models which need a spatial_file parameter
    elif model in ['SpatialMap','MapCubeFunction']:
        valid_pars=spatial_info.get('spatial_file') is not None

    #only valid model left has no required inputs
    else:
        valid_pars=True

    return valid_model and valid_pars

def valid_spectrum_input(spectrum_info):
    '''
    function to check if the input dict has valid parameter
    names and all necessary parameters for a Spectrum class object

    Parameters
    ----------
    spatial_info - dict
        dictionary with keys and values to create a
        Spectrum class object (see build_mode.model_components
        module for more inforamtion)

    Returns
    -------
    bool
        True if all parameter names are valid and all necessary
        parameters are specified, otherwise False
    '''
    
    #check if the model requested is valid
    model=spectrum_info.get('model')
    valid_model=model in spectral_models

    #check if necessary parameters are specified
    if model=='FileFunction':
        valid_pars=spectrum_info.get('spectrum_file') is not None

    #currently, all other models have default values
    #for all parameters
    else:
        valid_pars=True

    return valid_model and valid_pars
