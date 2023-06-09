#!/usr/bin/env python

from build_model.utilities import (
    angular_separation,
    get_ROI_from_event_file,
    build_region,
    valid_spatial_input,
    valid_spectrum_input)

from build_model.model_components import Spectrum,Spatial

import os

import warnings

import numpy as np

import pandas as pd

from functools import reduce

import astropy.io.fits as pyfits

from xml.dom import minidom

class SourceList:
    '''
    A class used to create a spatial and spectral XML source model for
    analysis of LAT data using the 4FGL catalog (data release versions 1, 2, or 3)

    ...

    Attributes
    ----------
    <initial>
    catalog_file : str
        path to 4FGL catalog file, FITS or XML format
    output_name : str
        path for output XML model file
    write_directory : str
        path for directory file will be saved in
    DR : int
        4FGL data release version of input catalog_file
    variability_threshold: float
        value to compare the variability index of a source when
        determining if it should be considered variable or not,
        depends on data release version
    ROI : list
        region of interest information for the model
        [right ascension, declination, radius]
    <subsequent>
    sources : dict
        dictionary of sources, keys are source names, selected
        from the 4FGL catalog which satisfy ROI and other requirements,
        will not exist until after the make_model method has
        been called
    num_extended : int
        number of extended sources added to the model, will not
        exist until after make_model method has been called and
        does not include diffuse emission components
    num_point : int
        number of point sources added to the model, will not
        exist until after make_model method has been called

    Methods
    -------
    add_point_source(source_name,RA,DEC,spectrum_model='PowerLaw',new_model_name=None,overwrite=False)
        adds a point source, not in 4FGL, at the given position
        to the XML model created via the make_model method
    add_source(source_name,spatial_info,spectrum_info,diffuse=False,new_model_name=None,overwrite=False)
        add a new source, not in 4FGL, with spatial and spectral information
        specified via input dictionaries, to an XML model created
        via the make_model method
    build_model()
        creates the spatial-spectral model, after sources have been
        selected from the 4FGL catalog, and saves the output XML file
    create_XML_model()
        calls methods to select sources from the 4FGL catalog to be included
        in the model and then calls a method to build the model
    get_sources_fits()
        selects sources from a FITS format version of the 4FGL catalog
    get_sources_xml()
        selects sources from a XML format version of the 4FGL catalog
    make_model(galactic_file="gll_iem_v07.fits",galactic_name='gll_iem_v07',
               isotropic_file="iso_P8R3_SOURCE_V3_v1.txt",isotropic_name='iso_P8R3_SOURCE_V3_v1',
               norms_free_only=False,extended_directory=None,free_radius=-1,max_free_radius=None,extra_radius=10,sigma_to_free=5,
               variable_free=True,force_point_sources=False,extended_catalog_names=False,make_region=True,region_file=None,
               galactic_index_free=True,use_old_names=False)
        checks inputs, calls methods to select sources and make the XML model,
        and makes a ds9-styel .reg file from the selected sources (if requested)
    Print()
        prints to screen details of the inputs (catalog, ROI)
        and outputs (file name and write directory)
    '''    

    def __init__(self,catalog_file,ROI,output_name='my_model.xml',DR=3,write_directory=''):
        '''
        Parameters
        ----------
        catalog: str
            path to FITS or XML version of 4FGL catalog
        ROI: list or str
            region of interest information as either a list [RA,DEC,radius]
            or the path to a Fermi LAT event file with the desired
            ROI selection in the header
        output_name: str
            desired file name of the output XML model, note that this
            should not be the full path
        DR: int
            data release version of the 4FGL catalog corresponding
            to the passed "catalog" parameter
        write_directory: str
            path to directory in which to save "output_name", if not specified
            will default to the current directory
        '''
        
        #some sanity checks on input values
        if not os.path.exists(catalog_file):
            raise FileNotFoundError(2,'Could not access catalog file',catalog_file)
        
        if os.path.exists(output_name):
            warnings.warn(f'Region XML model {output_name} already exists, will be overwritten\
 (note, some systems may be case-insensitive).')

        if DR not in [1,2,3,4]:
            raise ValueError(f'{DR = } is an invalid choice of data release, must be 1, 2, 3, or 4.')

        self.catalog_file=catalog_file
        self.output_name=os.path.join(write_directory,output_name)
        self.write_directory=write_directory
        self.DR=DR
        self.variability_threshold=27.69 if DR==4 \
                              else 24.725 if DR==3 \
                              else 21.666 if DR==2 \
                              else 18.475
        
        #self.fermi_dir=os.getenv('FERMI_DIR')

        #get the region of inerest information
        #either directly from a list passed in
        if isinstance(ROI,list):
            self.ROI=ROI

        #or from the header of an eventfile
        #passed in
        elif os.path.exists(ROI):
            self.ROI=get_ROI_from_event_file(ROI)

        #may need to rethink the way this else
        #statement is handled
        else:
            raise FileNotFoundError(2,'Could not access event file',ROI)
    
    
    def Print(self):
        '''
        print the inputs and desired output
        '''
        
        print('Catalog file: ',self.catalog_file)
        print('Output file name: ',self.output_name)
        print(f'Selecting {self.ROI[2]:.1f} degrees around (ra,dec)=({self.ROI[0]:.3f},{self.ROI[1]:.3f})')
    
    def make_model(self,galactic_file="gll_iem_v07.fits",galactic_name='gll_iem_v07',
               isotropic_file="iso_P8R3_SOURCE_V3_v1.txt",isotropic_name='iso_P8R3_SOURCE_V3_v1',
               norms_free_only=False,extended_directory=None,free_radius=-1,max_free_radius=None,extra_radius=10,sigma_to_free=5,
               variable_free=True,force_point_sources=False,extended_catalog_names=False,make_region=True,region_file=None,
               galactic_index_free=True,use_old_names=False):
        '''
        select sources from the 4FGL catalog (based on specified thresholds
        and logic), build and save an XML model, and optionally create a ds9-style
        .reg file
        all parameters are optional with reasonable defaults and all parameters
        become attributes, with the same name, of the SourceList object
        
        Parameters
        ----------
        galactic_file : str
            name of Galactic diffuse model file name, if full path is not
            specified will default to fermitools install location
        galactic_name : str
            name for Galactic diffuse component in XML model
        isotropic_file : str
            name of isotropic diffuse template file name, if full path is
            not specified will default to fermitools install lcoation
        isogtropic_name : str
            name for isotropic diffuse emission component in XML model
        norms_free_only : bool
            flag to only set spectral normalization parameters free to vary
            for any sources satisfying other requirements to be free
        extended_directory : str
            full path to directory with extended soure templates
            if not provided will use fermitools install location
        free_radius : float
            any sources further than this radius, in degrees, from ROI
            center will have all spectral parameters frozen to catalog
            values, unless they satisfy other requirements
        max_free_radius : float
            any sources further than this radius, in degrees, from ROI
            center will have all spectral parameters frozen to catalog
            values, regardless of other requirements
        extra_radius : float
            additional radius, in degrees, beyond ROI radius within which
            to include sources, with all parameters fixed, to account for
            the size of the point-spread function at low energies
        sigma_to_free : float
            sources with average significance values in the 4FGL catalog, if
            using the FITS version, will have all spectral parameters fixed
            at catalog values, though variability information may set the
            normalization parameters free.  if using the XML version of the
            catalog, this values corresponds to test statistic value
        variable_free : bool
            flag to set spectral normalization parameters free for sources
            which would otherwise be fixed if their variability index exceeds
            the corresponding DR version threshold, does not apply to sources
            outside the 'max_free_radius' boundary
        force_point_sources : bool
            flag to include extended sources, not including the Galactic and
            isotropic diffuse components, as point sources in the output XML model
        extended_catalog_names : bool
            flag to use the catalog names (e.g., 4FGL J####.#+####e) instead of
            association names for extended sources
        make_region : bool
            flag to make a ds9-style .reg file corresponding to the XML model,
            with color indicating free or fixed spectral parameter(s) and
            shape indicating spectral model
        region_file : str
            output name of ds9-style .reg file, if not supplied will construct
            to match name of output XML model file, this is assumed to be in the
            "write_directory", do not give the full path
        galactic_index_free : bool
            flag to modify the spectrum of the Galactic diffuse emission by a power law
            with free Index parameter
        use_old_names : bool
            flag to use the old naming convention where spaces are removed and
            an underscore is prepended to the source name (e.g., 4FGL J####.#+####
            becomes _4FGLJ####.#+####)
        '''

        if os.path.dirname(galactic_file)=='':
            #if self.fermi_dir is None:
            galactic_file=os.sep.join(['$(FERMI_DIR)','refdata','fermi','galdiffuse',galactic_file])
            #not sure if we should hardcode the path or not, makes it difficult if trying to work
            #with different fermitools conda installs
            #else:
            #    galactic_file=os.path.join(self.fermi_dir,'refdata','fermi','galdiffuse',galactic_file)

        if os.path.dirname(isotropic_file)=='':
            #if self.fermi_dir is None:
            isotropic_file=os.sep.join(['$(FERMI_DIR)','refdata','fermi','galdiffuse',isotropic_file])
            #not sure if we should hardcode the path or not, makes it difficult if trying to work
            #with different fermitools conda installs
            #else:
            #    isotropic_file=os.path.join(self.fermi_dir,'refdata','fermi','galdiffuse',isotropic_file)

        #create attributes associated with diffuse
        #emission components
        self.galactic_file=galactic_file
        self.galactic_name=galactic_name

        self.isotropic_file=isotropic_file
        self.isotropic_name=isotropic_name

        #create attributes associated with radial boundaries
        #and check values
        self.free_radius=(self.roi[2] if free_radius<=0 else free_radius)
        self.max_free_radius=(self.free_radius if max_free_radius is None else max_free_radius)
        
        if self.max_free_radius<self.free_radius:
            warnings.warn(f'max_free_radius ({max_free_radius:.1f}) is < free_radius ({free_radius:.1f}), making the parameter useless.')

        self.extra_radius=extra_radius

        #determine the extended source template directory to use
        self.extended_directory=extended_directory if extended_directory is not None else\
                    os.sep.join(['$(FERMI_DIR)','data','pyBurstAnalysisGUI','templates'])
        
        #create variables associated with boolean flags
        #less to pass in to methods this way
        #and can easily verify what was used
        self.variable_free=variable_free
        self.force_point_sources=force_point_sources
        self.extended_catalog_names=extended_catalog_names
        self.norms_free_only=norms_free_only
        self.use_old_names=use_old_names
        self.sigma_to_free=sigma_to_free
        self.make_region=make_region
        self.galactic_index_free=galactic_index_free

        #make the XML model
        self.create_XML_model()

        #make the ds9-style .reg model, if requested
        if make_region:
            #if name not specified, construct from output_name value
            if region_file is None or region_file=='':
                region_file=os.path.basename(self.output_name).split('.')[:-1]
                #allow that '.' might be in the name beyond just giving the extenstion
                region_file=reduce(lambda s1,s2:s1+s2,region_file).join(['ROI_','.reg'])

            #construct the full path to the file
            self.region_file=os.path.join(self.write_directory,region_file)

            #and create the .reg file
            print('Building ds9-style region file',end='...')
            build_region(self.region_file,self.sources)
            print(f'done!\nFile saved as {self.region_file}.')
    
    def create_XML_model(self):
        '''
        based on the value of output-name, calls either the get_sources_xml
        or get_sources fits method to select sources from the 4FGL catalog
        matching the ROI and other constraints and then calls the build_model
        method to create the XML model file
        '''
        
        print(f'Creating spatial and spectral model from the 4FGL DR-{self.DR} catalog: {self.catalog_file}.')
        
        if os.path.basename(self.catalog_file).split('.')[-1].lower()=='xml':
            self.get_sources_xml()
        else:
            self.get_sources_fits()

        #from either method we get a self.sources attribute which is a nested dictionary with all the sources
        self.build_model()

    def build_model(self):
        '''
        function to construct a spatial-spectral XML region model from
        the sources attribute of the SourceList object (requires get_sources_fits
        or get_sources_xml method to have been called first) also adds the
        isotropic and Galactic diffuse emission components
        '''
        
        #keep track of point and extended sources
        self.num_extended=0
        self.num_point=0

        #use the self.sources nested dictionary to create the output xml file
        output_xml=minidom.getDOMImplementation().createDocument(None,'source_library',None)
        output_xml.documentElement.setAttribute('title','source library')

        for source_name in self.sources.keys():
            #create a new 'source' element
            source_out=output_xml.createElement('source')

            #set basic attributes of the element
            source_out.setAttribute('name','_'+reduce(lambda s1,s2:s1+s2,source_name.split(' ')) if self.use_old_names else source_name)
            source_out.setAttribute('ROI_Center_Distance',f"{self.sources[source_name]['roi_distance']:.2f}")
            source_out.setAttribute('type','DiffuseSource' if self.sources[source_name]['Extended'] else 'PointSource')

            
            if self.sources[source_name]['Extended']:
                self.num_extended+=1
            else:
                self.num_point+=1
            
            #need to determine if parameters should be free or fixed
            #first, basic checks for conditions meaning everything is fixed
            if self.sources[source_name]['roi_distance']>self.max_free_radius or\
               (self.sources[source_name]['roi_distance']>self.free_radius and not self.sources[source_name]['variable']) or\
               (self.sources[source_name]['roi_distance']<=self.free_radius and not self.sources[source_name]['significant']\
                and not self.sources[source_name]['variable']):

                update_info=[]

                for spectrum_par in self.sources[source_name]['spectrum'].keys():
                    if spectrum_par not in ['model','spectrum_file']:
                        update_info.append((spectrum_par+'_free',False))

                self.sources[source_name]['spectrum'].update(update_info)
                del update_info


                self.sources[source_name].update([('free',False)])

            #now, check if the source is far from the ROI center but variable
            #or if close by but not formally significant or close by and significant but norms_only_free set
            #free the normalization, only
            elif (self.sources[source_name]['roi_distance']<=self.max_free_radius and\
                 self.sources[source_name]['roi_distance']>self.free_radius and self.sources[source_name]['variable']) or \
                 (self.sources[source_name]['roi_distance']<=self.free_radius and not self.sources[source_name]['significant']\
                  and self.sources[source_name]['variable']) or (self.sources[source_name]['roi_distance']<=self.free_radius\
                and self.sources[source_name]['significant'] and self.norms_free_only):

                update_info=[]
                
                for spectrum_par in self.sources[source_name]['spectrum'].keys():
                    if spectrum_par  not in ['model','spectrum_file']:
                        if spectrum_par in ['Prefactor','Integral','norm','Normalization']:
                            update_info.append((spectrum_par+'_free',True))
                        else:
                            update_info.append((spectrum_par+'_free',False))

                self.sources[source_name]['spectrum'].update(update_info)
                del update_info

                self.sources[source_name].update([('free',True)])

            #now, check if the source is close to the ROI center and significant
            #need to also consider if norms_free_only flag is set
            #this should be the only possiblity left, but for now we'll use elif and put an else
            #statement to catch things we haven't thought of, remove after successful testing
            elif self.sources[source_name]['roi_distance']<=self.free_radius and self.sources[source_name]['significant']:

                update_info=[]
                
                for spectrum_par in self.sources[source_name]['spectrum'].keys():
                    if spectrum_par not in ['model','spectrum_file']:
                        if spectrum_par in ['Scale','Index2','Eb']:
                            update_info.append((spectrum_par+'_free',False))
                        else:
                            update_info.append((spectrum_par+'_free',True))

                self.sources[source_name]['spectrum'].update(update_info)
                del update_info

                self.sources[source_name].update([('free',True)])
                
            #the following else statement is for testing only
            #should not drop into that section of code, remove after many successful tests
            #replace the 'elif' above with simply 'else'
            else:
                print('Should not have arrived here! Look at nested dictionary to see what happened.')
                print(self.sources[source_name])
                raise RuntimeError()

            #now we need to add the extended_directory to extended source file names
            if self.sources[source_name]['Extended'] and\
               self.sources[source_name]['spatial']['spatial_model'] in ['SpatialMap','MapCubeFunction']:
                self.sources[source_name]['spatial']['spatial_file']=os.path.join(self.extended_directory,
                                        self.sources[source_name]['spatial']['spatial_file'])
                print(f'Extended source {source_name} in model, make sure that the following path to extended template is correct:',end=' ')
                print(self.sources[source_name]['spatial']['spatial_file'])

            #now that we have added flags for spectral parameters to be free or fixed
            #and possibly adjusted extended source file paths
            #let's create the Spectrum and Spatial objects (which will add
            #spectrum and spatialModel elements to the source element
            #and then add them to the source and add that to the larger XML document
            spectrum=Spectrum(**self.sources[source_name]['spectrum'])
            source_out.appendChild(spectrum.spectrum)

            spatial=Spatial(**self.sources[source_name]['spatial'])
            source_out.appendChild(spatial.spatial)

            #now add the source to the output_xml
            output_xml.documentElement.appendChild(source_out)

        #we've cycled through all of the sources, now we need to add the diffuse components
        #first the Galactic diffuse
        Galactic=output_xml.createElement('source')
        Galactic.setAttribute('name',self.galactic_name)
        Galactic.setAttribute('type','DiffuseSource')

        galactic_spectrum=Spectrum(model='PowerLaw',Prefactor=1,Index=0,Scale=100,Prefactor_free=True,
                       Index_free=self.galactic_index_free)
        galactic_spatial=Spatial(spatial_model='MapCubeFunction',spatial_file=self.galactic_file)

        Galactic.appendChild(galactic_spectrum.spectrum)
        Galactic.appendChild(galactic_spatial.spatial)

        output_xml.documentElement.appendChild(Galactic)

        #now the isotropic component
        Isotropic=output_xml.createElement('source')
        Isotropic.setAttribute('name',self.isotropic_name)
        Isotropic.setAttribute('type','DiffuseSource')

        isotropic_spectrum=Spectrum(model='FileFunction',spectrum_file=self.isotropic_file,apply_edisp='false')
        isotropic_spatial=Spatial(spatial_model='ConstantValue')

        Isotropic.appendChild(isotropic_spectrum.spectrum)
        Isotropic.appendChild(isotropic_spatial.spatial)
        
        output_xml.documentElement.appendChild(Isotropic)

        #put all of this in a string and write to output file
        out_string=filter(lambda s: len(s) and not s.isspace(),
                  output_xml.toprettyxml(' ').splitlines(True))

        with open(self.output_name,'w') as output_file:
            output_file.write(''.join(out_string))

        #print information about the model for the user
        if self.force_point_sources:
            print(f'Added {self.num_point} point sources, note that any extended sources\
 (except diffuse components) have been modeled as point sources.')

        else:
            print(f'Added {self.num_point} point sources and {self.num_extended} extended sources.')

            if self.num_extended>0:
                print('If you plan to analyze LAT data with unbinned likelihood,\
 you will need to run gtdiffrsp for the extended sources or rerun makeModel\
 with force_point_sources set to True.')


    def get_sources_xml(self):
        '''
        selects sources from a XML version of the 4FGL catalog
        which satisfy the ROI and other requirements
        '''

        #read in the catalog
        catalog=minidom.parse(self.catalog_file)
        catalog_sources=np.array(catalog.getElementsByTagName('source'))

        #traverse through the sources once, just to get positional information
        right_ascensions,declinations=[],[]
        for source in catalog_sources:
            #try to get information potentially from the 'header'
            try:
                right_ascensions.append(float(source.getAttribute('RA')))
                declinations.append(float(source.getAttribute('DEC')))

            #if that doesn't work, access the spatialModel element directly
            except:
                for parameter in source.getElementsByTagName('spatialModel')[0].getElementsByTagName('parameter'):
                    if parameter.getAttribute('name')=='RA':
                        right_ascensions.append(float(parameter.getAttribute('value')))

                    elif parameter.getAttribute('name')=='DEC':
                        declinations.append(float(parameter.getAttribute('value')))

        #recast the lists as numpy arrays for use in
        #angular separation function
        right_ascensions=np.array(right_ascensions)
        declinations=np.array(declinations)

        #get distance of each source from the ROI center
        self.source_distances=angular_separation(self.ROI[0],self.ROI[1],
                             right_ascensions,declinations)

        #now we want to save time by only getting sources
        #satisfying the ROI requirements
        source_mask=self.source_distances<=(self.ROI[2]+self.extra_radius)
        self.source_distances=self.source_distances[source_mask]

        #sort the distances in increasing order
        distance_index=self.source_distances.argsort()
        self.source_distances=self.source_distances[distance_index]

        #apply the masking and sorting to the catalog sources
        #right ascensions, and declinations
        catalog_sources=catalog_sources[source_mask][distance_index]
        right_ascensions=right_ascensions[source_mask][distance_index]
        declinations=declinations[source_mask][distance_index]

        #now traverse the shorter list of sources to create the self.sources nested dictionary
        self.sources={}
        
        for source,RA,DEC,distance in zip(catalog_sources,right_ascensions,declinations,self.source_distances):
            #get the name and position info
            name=source.getAttribute('name')
            self.sources.update([(name,{'roi_distance':distance,'RA':RA,'DEC':DEC})])

            #check if the source is extended or point-like
            if source.getAttribute('type')=='DiffuseSource':
                self.sources[name].update([('Extended',not self.force_point_sources)])
            else:
                self.sources[name].update([('Extended',False)])

            #evaluate variability and significance threshold info
            self.sources[name].update([('variable',float(source.getAttribute('Variability_Index'))>=self.variability_threshold),
                ('significant',float(source.getAttribute('TS_value'))>=self.sigma_to_free)])

            #now get the spectrum information
            spectrum=source.getElementsByTagName('spectrum')
            
            self.sources[name].update([('spectrum',
                 {'model':spectrum[0].getAttribute('type')})])
            
            for parameter in spectrum[0].getElementsByTagName('parameter'):
                self.sources[name]['spectrum'].update([(parameter.getAttribute('name'),
                float(parameter.getAttribute('value'))*float(parameter.getAttribute('scale')))])
            
            #now get the spatial information
            spatial=source.getElementsByTagName('spatialModel')

            #more involved for extended sources
            if self.sources[name]['Extended']:
                self.sources[name].update([('spatial',
                    {'spatial_model':spatial[0].getAttribute('type')})])

                if self.sources[name]['spatial'].get('spatial_model') in ['SpatialMap','MapCubeFunction']:
                    self.sources[name]['spatial'].update([('spatial_file',
                                        spatial[0].getAttribute('file').split('/')[-1])])

                for parameter in spatial[0].getElementsByTagName('parameter'):
                    self.sources[name]['spatial'].update([(parameter.getAttribute('name'),
                    parameter.getAttribute('value'))])

            #easier for point sources
            else:
                self.sources[name].update([('spatial',
                    {'spatial_model':'SkyDir','RA':RA,'DEC':DEC})])
            
            

    def get_sources_fits(self):
        '''
        selects sources from a FITS version of the 4FGL catalog
        which satisfy the ROI and other requirements
        '''
        
        #first, open the catalog FITS file and extract the information
        #we care about, putting the info into a couple of pandas dataframes
        #for ease of dealing with and searching through
        with pyfits.open(self.catalog_file) as catalog:
            #we're interested in only two of the extensions
            source_info=catalog['LAT_Point_Source_Catalog'].data.field
            extended_info=catalog['ExtendedSources'].data.field

            #access the necessary info for extended sources
            #stripping extraneous white spaces when necessary
            extended_sources=pd.DataFrame(np.c_[[name.strip() for name in extended_info('Spatial_Filename')],
                                                [function.strip() for function in extended_info('Spatial_Function')],
                                                extended_info('Model_SemiMajor'),
                                                extended_info('RAJ2000'),
                                                extended_info('DEJ2000')],
                                          columns=['file','function','extent','RA','DEC'],
                                          index=[name.strip() for name in extended_info('Source_Name')])

            #get the necessary information from the main catalog extension
            #stripping extraneous white spaces when necessary
            #and accounting for differences in DR versions
            catalog_sources=pd.DataFrame(np.c_[source_info('RAJ2000'),
                                               source_info('DEJ2000'),
                                               source_info('Signif_Avg'),
                                               source_info('Variability_Index'),
                                               [name.strip() for name in source_info('Extended_Source_Name')],
                                               source_info('Pivot_Energy'),
                                               source_info('PL_Flux_Density'),
                                               source_info('PL_Index'),
                                               source_info('LP_Flux_Density'),
                                               source_info('LP_Index'),
                                               source_info('LP_beta'),
                                               source_info('PLEC_Flux_Density'),
                                               source_info(('PLEC_IndexS' if self.DR==3 else 'PLEC_Index')),
                                               source_info(('PLEC_ExpfactorS' if self.DR==3 else 'PLEC_Expfactor')),
                                               source_info('PLEC_Exp_Index'),
                                               [model.strip() for model in source_info('SpectrumType')]],
                                         columns=['RA','DEC','Signif_Avg','Variability_Index',
                                                  'Extended_Name','Pivot_Energy','PL_Flux_Density',
                                                  'PL_Index','LP_Flux_Density','LP_Index',
                                                  'LP_beta','PLEC_Flux_Density','PLEC_Index1',
                                                  'PLEC_Expfactor','PLEC_Index2','Model'],
                                         index=[name.strip() for name in source_info('Source_Name')])

        #now we need to convert some columns to floats
        to_float_columns=['RA','DEC','Signif_Avg','Variability_Index',
                          'Pivot_Energy','PL_Flux_Density','PL_Index',
                          'LP_Flux_Density','LP_Index','LP_beta',
                          'PLEC_Flux_Density','PLEC_Index1',
                          'PLEC_Expfactor','PLEC_Index2']
        
        catalog_sources[to_float_columns]=catalog_sources[to_float_columns].astype(float)
        
        #calculate the distance from the roi for each source
        catalog_sources=pd.concat([catalog_sources,
                                   pd.DataFrame(angular_separation(self.ROI[0],self.ROI[1],
                                                    catalog_sources.RA.to_numpy(),
                                                    catalog_sources.DEC.to_numpy()),columns=['roi_distance'],
                                                index=catalog_sources.index)],axis=1).copy()

        #trim things down to just those sources we want to add
        catalog_sources=catalog_sources[catalog_sources.roi_distance<=(self.ROI[2]+self.extra_radius)]

        #and now sort the columns in order of ascending roi_distance
        catalog_sources.sort_values(by='roi_distance',ascending=True,inplace=True)

        #save the distances and make the sources dictionary
        self.source_distances=catalog_sources.roi_distance.to_numpy(copy=True)

        self.sources={}

        #cycle through the sources
        for source_name,row in catalog_sources.iterrows():
            #check if the source is extended or not, decide on name for source
            #and if it is extended
            if source_name[-1]=='e':
                name=source_name if self.extended_catalog_names else row.Extended_Name
                self.sources.update([(name,{'roi_distance':row.roi_distance})])
                self.sources[name].update([('Extended',not self.force_point_sources)])

            else:
                name=source_name
                self.sources.update([(name,{'roi_distance':row.roi_distance})])
                self.sources[name].update([('Extended',False)])
            
            #evaluate variability and significance threshold info
            #as well as the RA and DEC to be used if building a .reg file
            self.sources[name].update([('variable',row.Variability_Index>=self.variability_threshold),
                         ('significant',row.Signif_Avg>=self.sigma_to_free),
                         ('RA',row.RA),('DEC',row.DEC)])

            #now get source spectral information in the form we need
            #I need to verify that PowerLaw, LogParabola, and PLSuperExpCutoff4
            #are the only models we have in 4FGL DR 1, 2, and 3
            if row.Model=='PowerLaw':
                self.sources[name].update([('spectrum',
                    {'model':'PowerLaw',
                     'Prefactor':row.PL_Flux_Density,
                     'Index':row.PL_Index,
                     'Scale':row.Pivot_Energy})])

            elif row.Model=='LogParabola':
                self.sources[name].update([('spectrum',
                    {'model':'LogParabola',
                     'norm':row.LP_Flux_Density,
                     'alpha':row.LP_Index,
                     'beta':row.LP_beta,
                     'Eb':row.Pivot_Energy})])

            else:
                if self.DR==3:
                    self.sources[name].update([('spectrum',
                        {'model':'PLSuperExpCutoff4',
                         'Prefactor':row.PLEC_Flux_Density,
                         'IndexS':row.PLEC_Index1,
                         'Scale':row.Pivot_Energy,
                         'ExpfactorS':row.PLEC_Expfactor,
                         'Index2':row.PLEC_Index2})])

                else:
                    #for PLSuperExpCutoff2, we need to modify the flux density value in the catalog
                    #do it here and not in the PLSuperExpCutoff2 Spectrum class function
                    self.sources[name].update([('spectrum',
                        {'model':'PLSuperExpCutoff2',
                         'Prefactor':row.PLEC_Flux_Density*np.exp(row.Expfactor*row.Pivot_Energy**row.PLEC_Index2),
                         'Index1':row.PLEC_Index1,
                         'Scale':row.Pivot_energy,
                         'Expfactor':row.PLEC_Expfactor,
                         'Index2':row.PLEC_Index2})])

            #now for the spatial information
            #check if extended, this flag already uses force_point_sources flag
            if self.sources[name]['Extended']:
                #do extended stuff
                ext_row=extended_sources.loc[row.Extended_Name]

                #if one of the 'Radial' models
                if ext_row.function[:6]=='Radial':
                    self.sources[name].update([('spatial',
                        {'spatial_model':'RadialGaussian' if ext_row.function=='RadialGauss' else ext_row.function,
                         'RA':ext_row.RA,
                         'DEC':ext_row.DEC})])

                    if ext_row.function=='RadialDisk':
                        self.sources[name]['spatial'].update([('Radius',ext_row.extent)])

                    else:
                        self.sources[name]['spatial'].update([('Sigma',float(ext_row.extent)/(-2*np.log(0.32))**0.5)])

                #otherwise it is a spatial map
                else:
                    self.sources[name].update([('spatial',
                        {'spatial_model':'SpatialMap',
                         'spatial_file':ext_row.file})])

            #point source
            else:
                self.sources[name].update([('spatial',
                    {'spatial_model':'SkyDir',
                     'RA':row.RA,
                     'DEC':row.DEC})])

    def add_source(self,source_name,spatial_info,spectrum_info,
                   diffuse=False,new_model_name=None,overwrite=False):
        '''
        method to add a source, not in the 4FGL catalog, to an existing
        XML model

        Parameters
        ----------
        source_name : str
            name to be given to source which will be added to
            the XML model
        spatial_info : dict
            dictionary with necessary info to create a Spatial class object,
            see build_model.model_components module for more details
        spectrum_info : dict
            dictionary with necessary info to create a Spectrum class object,
            see build_model.model_components module for more details
        diffuse : bool
            flag to indicate if the source is extended (True) or
            point-like (False)
        new_model_name : str
            file name for new XML model, if not specified the existing
            model will be overwritten if the overwrite parameter is True
        overwrite : bool
            flag to overwrite existing file of same name
        '''
        
        #first, do a check on the model to make sure it exists
        if not os.path.exists(self.output_name):
            raise RuntimeError(f'{self.output_name} has not been created yet,\
 cannot use add_source until make_model has been successfully run.')

        #assume the model is in the same directory as the
        #existing model, add the full path information
        if new_model_name is not None:
            if os.path.dirname(new_model_name)=='':
                new_model_name=os.path.join(\
                    os.path.dirname(self.output_name),
                    new_model_name)
                warnings.warn('Assuming same save directory as main model.')

            if os.path.exists(new_model_name):
                if overwrite:
                    warnings.warn(f'File {new_model_name} exists, will be overwritten.')
                else:
                    raise RuntimeError(f'File {new_model_name} exists but overwrite flag set to False.')

        #if a name isn't provided, we will overwrite the existing
        #file, but only if overwrite is True
        else:
            if overwrite:
                new_model_name=self.output_name
                warnings.warn('No new file name given, will overwrite existing model.')

            else:
                raise RuntimeError('No new file name given but overwrite flag set to False.')

        #need to ensure that the dictionaries include minimum info
        #before proceeding with trying to add a source
        #note, these functions do not check if extra parameters
        #are included, which will cause an error
        if not valid_spatial_input(spatial_info):
            raise RuntimeError('Dictionary with spatial info missing required\
 inputs or has invalid model selection.')

        if not valid_spectrum_input(spectrum_info):
            raise RuntimeError('Dictionary with spectrum info missing required\
 inputs or has invalid model selection.')

        #make the Spatial and Spectrum objects which
        #create the necessary XML elements
        spatial=Spatial(**spatial_info)
        spectrum=Spectrum(**spectrum_info)

        #now create a new xml
        output_xml=minidom.getDOMImplementation().createDocument(None,'source_library',None)
        output_xml.documentElement.setAttribute('title','source library')

        #make a source element
        source_out=output_xml.createElement('source')
        source_out.setAttribute('name',source_name)
        
        #get the angular separation if spatial model has RA and DEC info
        if {'RA','DEC'}.issubset(spatial_info.keys()):
            source_sep=angular_separation(self.ROI[0],self.ROI[1],
                                spatial_info.get('RA'),spatial_info.get('DEC'))
            source_out.setAttribute('ROI_Center_Distance',f'{source_sep:.2f}')
        
        source_out.setAttribute('type','DiffuseSource' if diffuse else 'PointSource')

        #append the spectral and spatial elements
        source_out.appendChild(spectrum.spectrum)
        source_out.appendChild(spatial.spatial)

        #add the source element to our new xml
        output_xml.documentElement.appendChild(source_out)

        #now cycle through sources already in the model and
        #add them to the new XML document
        current_model=minidom.parse(self.output_name)
        current_sources=np.array(current_model.getElementsByTagName('source'))

        for source in current_sources:
            output_xml.documentElement.appendChild(source)

        #put it all in a string and write out
        out_string=filter(lambda s: len(s) and not s.isspace(),
                  output_xml.toprettyxml(' ').splitlines(True))

        with open(new_model_name,'w') as output_file:
            output_file.write(''.join(out_string))

        print(f'{source_name} added to model, saved as {new_model_name}')

        self.output_name=new_model_name
        print('SourceList object attribute "output_name" now points to new file')

    #a convenience method to add a point source easily
    #will call the add_source method
    #spectrum_model can be the name of the model or a dictionary with
    #the necessary model info
    def add_point_source(self,source_name,RA,DEC,spectrum_model='PowerLaw',new_model_name=None,overwrite=False):
        '''
        method to add a point source, not in the 4FGL catalog,
        to an existing XML model, this is a convenience method
        which will, itself, call the add_source method

        Parameters
        ---------
        source_name : str
            name to be given to source which will be added to
            the XML model
        spectrum_model : str or dict
            either the name of the spectral model for the source or a
            dictionary with necessary info to create a Spectrum class object,
            see build_model.model_components module for more details
        new_model_name : str
            file name for new XML model, if not specified the existing
            model will be overwritten if the overwrite parameter is True
        overwrite : bool
            flag to overwrite existing file of same name
        '''

        #determine what type of spectrum_model input we have
        if isinstance(spectrum_model,dict):
            spectrum_info=spectrum_model

        else:
            spectrum_info={'model':spectrum_model}

        #call the add_source method
        self.add_source(source_name,
                        {'RA':RA,'DEC':DEC,'spatial_model':'SkyDir'},
                         spectrum_info,diffuse=False,
                         new_model_name=new_model_name,overwrite=overwrite)

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
        the corresponding True or False to Input
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
    
if __name__=='__main__': cli()
