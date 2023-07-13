#!/usr/bin/env python

#this version of the script is designed to take advantage
#of the planned, but not yet implemented, lat_source_model
#package installable via pip
#previous interactive python session code would import from
#this package and not make4FGLxml, this would remain
#as a command line interface
from LATSourceModel import SourceList

def mybool(Input):
    '''
    function to take command line input meant for boolean
    arguments and go from a set of possible strings
    to the correct value so that every input doesn't
    simply evaluate to True

    Parameters
    ----------
    Input - str
        input from command line, interpreted as a string

    Returns
    -------
    bool
        the True or False bool matching "Input"
    '''
    
    return {'True':True,'False':False,
            'T':True,'F':False,
            't':True,'f':False,
            'TRUE':True,'FALSE':False,
            "true":True,"false":False,
            "1":True,"0":False}.get(Input)

#function to be done when called from command line interface
def cli():
    '''
    function to take command line input, create a SourceList object
    and make a XML model using supplied ROI and requirements.  See
    the helpString below and associated argument help strings for more
    information (e.g., make4FGLxml.py --help)
    '''
    
    import argparse
    
    helpString="Creates an XML model from the 4FGL catalog (DR 1, 2, or 3) using FITS or XML catalog version\
            for a specific ROI specified through input coordinates or taken from an input event file,\
            the default radius for including sources is 10 degrees beyond the ROI radius,\
            sources with free parameters within the original extraction radius are chosen\
            based on nearness to center, significance, and variability."
    parser=argparse.ArgumentParser(description=helpString)
    
    parser.add_argument("catalog",type=str,help="Path to catalog file to use, can be FITS or xml.")
    
    parser.add_argument("-ev","--event_file",type=str,default=None,help="Event file with ROI information in header.")
    
    parser.add_argument("--RA",type=float,default=None,
                help="Right Ascension (J2000) in degrees of ROI center, only used if no event file given.")
    
    parser.add_argument("--DEC",type=float,default=None,
                help="Declination (J2000) in degrees of ROI center, only used if no event file given.")
    
    parser.add_argument("--radius",type=float,default=None,
                help="Radius, in degrees, of ROI, only used if no event file given.")
    
    parser.add_argument("-o","--output_name",type=str,default='my_model.xml',
                help="Name of output xml file, is set to overwrite files of same name.")
    
    parser.add_argument("-G","--galactic_file",type=str,default='gll_iem_v07.fits',
                help="Name and location of Galactic diffuse model to use, will default to gll_iem_v07.fits.")
    
    parser.add_argument("-g","--galactic_name",type=str,default='gll_iem_v07',
                help="Name of Galactic diffuse component in output model, will default to gll_iem_v07.")
    
    parser.add_argument("-I","--isotropic_file",type=str,default='iso_P8R3_SOURCE_V3_v1.txt',
                help="Name of isotropic diffuse template for output model, will default to P8R3 SOURCE class model.")
    
    parser.add_argument("-i","--isotropic_name",type=str,default='iso_P8R3_SOURCE_V3_v1'
                ,help="Name of isotropic diffuse component in output model, default is for P8R3 SOURCE class.")
    
    parser.add_argument("-N","--norms_free_only",type=mybool,default=False
                ,help="Flag to only let the normalizations of parameters be free, default is False.",
                nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
    
    parser.add_argument("-e","--extended_directory",type=str,default='',
                help="Path to directory with LAT extended source templates, will default to STs default.")
    
    parser.add_argument("-r","--free_radius",type=float,default=-1.,
                help="Radius, in degrees, from ROI center beyond which all source parameters should be fixed,\
 will default to selection radius.")
    
    parser.add_argument("-R","--max_free_radius",type=float,default=None,
                help="Absolute maximum radius, in degrees, from ROI center beyond which all source parameters should be fixed,\
 even variable sources will not be freed beyond this radius, defaults to free_radius value.")
    
    parser.add_argument("-ER","--extra_radius",type=float,default=10.,
                help="Radius beyond event file ROI out to which sources will be included in the model with all parameters fixed,\
 default is 10, good for analyses starting around a few hundred MeV, can be decreased for high energy only fits.")
    
    parser.add_argument("-s","--sigma_to_free",type=float,default=5.,
                help="Average significance below which all source parameters are fixed, defaults to 5.\
 Note, if using the 3FGL catalog xml file as input, this is actually a cut on TS, so adjust accordingly.")
    
    parser.add_argument("-v","--variable_free",type=mybool,default=True,
                help="Flag to set normalization of significantly variable sources,\
 even if source is beyond radius limit or below TS limit, default is True.",
                choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
    
    parser.add_argument("-p","--force_point_sources",type=mybool,default=False,
                help="Flag to cast extended sources as point sources, default is False.",
                nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
    
    parser.add_argument("-E2C","--extended_catalog_names",type=mybool,default=False,
                help="Flag to use catalog names for extended sources, default is False.",
                nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
    
    parser.add_argument("-m","--make_region",type=mybool,default=True,
                help="Flag to create ds9 region file as well as the xml model, default is True.",
                choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])

    parser.add_argument('-rf','--region_file',type=str,default='',
                help='Name of output .reg file, if make_region set to true.  Will default to\
 ROI_"output_name".reg if not specified.')
    
    parser.add_argument("-GIF","--galactic_index_free",type=mybool,default=True,
                help="Flag to use a power-law modification to the Galactic diffuse model spectrum\
 and have the index be free, default is True.",
                nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
    
    parser.add_argument("-wd","--write_directory",type=str,default='',
                help="Directory to write the output ds9 region file in if not the current working directory\
 or if you are specifying the full path to the newly made XML file.")
    
    parser.add_argument("-ON","--use_old_names",type=mybool,default=False,
                help="Flag to use the make2FLGxml style naming convention, underscore before name and no spaces,\
 default is False.",nargs="?",const=True,
                choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
    
    parser.add_argument('-DR',type=int,default=3,help='Choice of data release, 3 for 4FGL-DR3, 2 for 4FGL-DR2, and 1 for 4FGL.',
                choices=[1,2,3])

    parser.add_argument('--new',type=str,default=None,help='RA, DEC, source_name, and model for a non-4FGL point source to be\
 added to the model.  Must be a string contained within single quotes with keys in double quotes.')
    
    args=parser.parse_args()

    #determine how the ROI information will be derived,
    #event file or location and radius specified
    if args.event_file is None:
        if args.RA is None or args.DEC is None or args.radius is None:
            raise ValueError('No Fermi LAT event file provided but one of RA, DEC, or radius is None.')
        ROI=[args.RA,args.DEC,args.radius]
    else:
        ROI=args.event_file

    #Make a SourceList object and then call the makeModel method
    source_list=SourceList(args.catalog,ROI,args.output_name,args.DR,args.write_directory)
    
    source_list.make_model(args.galactic_file,args.galactic_name,
                  args.isotropic_file,args.isotropic_name,
                  args.norms_free_only,args.extended_directory,
                  args.free_radius,args.max_free_radius,args.extra_radius,
                  args.sigma_to_free,args.variable_free,args.force_point_sources,
                  args.extended_catalog_names,args.make_region,args.region_file,
                  args.galactic_index_free,args.use_old_names)

    #check if a non-4FGL source should be added
    if args.new is not None:
        import json

        #try to load the input dictionary, be ready to catch an exception
        #if the correct usage of ' and " is not followed and fix it
        try:
            new_source=json.loads(args.new)

        #swap out single quotes for double quotes
        except json.decoder.JSONDecodeError:
            alter_new=args.new.replace("'",'"')
            new_source=json.loads(alter_new)

        #make sure that we have the 4 required keys
        if {'RA','DEC','source_name','model'}.issubset(new_source.keys()):

            #use the add_point_source method, assume that since things are
            #being called from the script we should overwrite the XML model
            #and .reg files just made with the versions including the new source
            source_list.add_point_source(source_name=new_source.get('source_name'),
                                         RA=new_source.get('RA'),
                                         DEC=new_source.get('DEC'),
                                         spectrum_model=new_source.get('model'),
                                         update_reg=args.make_region,
                                         overwrite=True)

        else:
            print(f'Information for non-4FGL point source provided does not have the\
 correct keys\nExpected: ["RA","DEC","source_name","model"]\nReceived: {list(new_source.keys())}.')
            print('Cannot add new point source.')

    
if __name__=='__main__': cli()
