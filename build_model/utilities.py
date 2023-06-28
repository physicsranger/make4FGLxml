import astropy.io.fits as pyfits
import numpy as np

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

#function to extract the region of interest information from
#a Fermi LAT event FITS file
def get_ROI_from_event_file(event_file_name):
    #open the file and get the DSKEY info from the header
    with pyfits.open(event_file_name) as event_file:
        num=event_file[1].header['NDSKEYS']
        header=event_file[1].header
        
        dstypes={event_file[1].header[f'DSTYP{idx+1}']:\
            event_file[1].header[f'DSVAL{idx+1}']\
            for idx in range(event_file[1].header['NDSKEYS'])}
    
    #now attempt to get the position information and return
    try:
        RA,DEC,radius=dstypes.pop('POS(RA,DEC)').lower().strip('circle()').split(',')
        return float(RA),float(DEC),float(radius)
    
    except KeyError:
        print(f"No position keyword found in header of file {event_file_name} (assuming position is RA and DEC).")
        return None,None,None


#function to calculate the angular separation between a reference location
#on the sky and either single location on the sky or an array of sky locations
#takes arguments in degrees, returns values in degrees
def angular_separation(ref_RA,ref_DEC,RA_values,DEC_values):
    #for ease, create a degree-to-radians conversion variable
    d2r=np.pi/180
    
    #round the number from the calculation to avoid issues with the
    #arccos function
    return np.arccos((np.cos(ref_DEC*d2r)\
        *np.cos(DEC_values*d2r)*np.cos((ref_RA-RA_values)*d2r)\
        +np.sin(ref_DEC*d2r)*np.sin(DEC_values*d2r)).round(6))/d2r
    
#function to create a ds9 region file
def build_region(region_file_name,sources):
    #build a shape dictionary and a color list
    shape={'PLSuperExpCutoff':'cross','LogParabola':'diamond',
        'PowerLaw':'circle'}
    color=('magenta','green')
    
    with open(region_file_name,'w') as region_file:
        region_file.write('# Region file format: DS9\n')
        region_file.write('# Created by make4FGLxml\n')
        region_file.write('global font="roman 10 normal" move =0\n')
        
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
            entry+=f"text={ {source_name} }\n"
            region_file.write(entry)

def valid_spatial_input(spatial_info):
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
    #check if the model requested is valid
    model=spectrum_info.get('model')
    valid_model=mode in spectral_models

    #check if necessary parameters are specified
    if model=='FileFunction':
        valid_pars=spectrum_info.get('spectrum_file') is not None

    #currently, all other models have default values
    #for all parameters
    else:
        valid_pars=True

    return valid_model and valid_pars
