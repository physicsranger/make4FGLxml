#!/usr/bin/env python

from build_model.utilities import (
        angular_separation,
        get_ROI_from_event_file,
        build_region)

from build_model.model_components import(
        PowerLaw,
        PLSuperExpCutoff2,
        PLSuperExpCutoff4,
        LogParabola,
        FileFunction,
        SkyDir,
        Radial,
        SpatialMap,
        ConstantValue,
        MapCubeFunction)

import os,warnings
import numpy as np
import pandas as pd
from functools import reduce
import astropy.io.fits as pyfits
from xml.dom import minidom

class SourceList:
	#arguments are:
	#catalog_file (string, filename of LAT source list fits file in catalog format)
	#ROI (either a list with [RA,DEC,radius] or string with name of event file to extract ROI info from)
	#output_name (string, name of output xml file, defaults to my_model.xml)
        #DR (int, data release version, default is 3)
        def __init__(self,catalog_file,ROI,output_name='my_model.xml',DR=3):
                #some sanity checks on input values
                if not os.path.exists(catalog_file):
                        raise FileNotFoundError(2,'Could not access catalog file',catalog_file)
                
		if os.path.exists(output_name):
                        warnings.warn(f'Region model {output_name} already exists, will be overwritten.')

                if DR not in [1,2,3]:
                        raise ValueError(f'DR={DR} is an invalid choice of data release, must be 1, 2, or 3.')

                self.catalog_file=catalog_file
                self.output_name=output_name
                self.DR=DR
                self.variability_threshold=24.725 if DR==3 else 21.666 if DR==1 else 18.475

                if isinstance(ROI,list):
                        self.ROI=ROI
                        
                elif os.path.exists(ROI)
                        self.ROI=get_ROI_from_event_file(ROI)
                else:
                        raise FileNotFoundError(2,'Could not access event file',ROI)
        
	#define a quick print function to make sure everything looks right
	def Print(self):
		print('Catalog file: ',self.catalog_file)
		print('Output file name: ',self.output_name)
		print(f'Selecting {self.roi[2]:.1f} degrees around (ra,dec)=({self.roi[0]:.3f},{self.roi[1]:.3f})')
	
	#make the xml file
	#arguments are:
	#galactic_file (str) -- optional, location and name of Galactic diffuse model to use
	#galactic_name (str) -- optional, name of Galactic diffuse component to use in xml model
	#isotropic_file (str) -- optional, location and name of Isotropic diffuse template to use
	#isotropic_name (str) -- optional, name of Isotropic diffuse component to use in xml model
	#norms_free_only (bool) -- optional, flag to only set normalizations parameters free
	#extended_directory (str) -- optional, directory with extended source templates
	#free_radius (float) -- optional, radius in degrees from center of ROI beyond which source parameters are fixed
	#max_free_radius (float) -- optional, absolute maximum radius beyond which sources are fixed, this may be necessary when doing binned analysis and a variable source beyond free_radius would be set free but this source is beyond the boundaries of the square region used for the binned likelihood
	#extra_radius (float) -- optional, radius beyond ROI radius in event file out to which sources will be included with fixed parameters, defaul tof 10 is good for analyses starting around 100 MeV, but for higher energy fits this can be decreased
	#sigma_to_free (float) -- optional, average significance, using FITS catalog file, below which source parameters are fixed, even if within free_radius.  This corresponds to TS_value if using the XML catalog file.
	#variable_free (float) -- optional, variability index above which source parameters are free, if beyond free_radius and/or below sigma_to_free only the normalization parameter is set free
	#force_point_sources (bool) -- optional, flag to force extended sources to be point sources
	#extended_catalog_names (bool) -- optional, flag to force use catalog names for extended sources (only matters if using catalog FITS file)
	#make_region (bool) -- optional, flag to also generate ds9 region file
	#galactic_index_free (bool) -- optional, the Galactic diffuse is given a power-law spectral shape but the by default the index is frozen, setting this flag to True allows that to be free for additional freedom in diffuse fit
	#use_old_names (bool) -- optional, flag to use the naming convention from make1FGLxml.py and make2FGLxml.py with a leading underscore and no spaces
	def make_model(self,galactic_file="$(FERMI_DIR)/refdata/fermi/galdiffuse/gll_iem_v07.fits",galactic_name='gll_iem_v07',
                      isotropic_file="$(FERMI_DIR)/refdata/fermi/galdiffuse/iso_P8R3_SOURCE_V3_v1.txt",isotropic_name='iso_P8R3_SOURCE_V3_v1',
                      norms_free_only=False,extended_directory=None,free_radius=-1,max_free_radius=None,extra_radius=10,sigma_to_free=5,
                      variable_free=True,force_point_sources=False,extended_catalog_names=False,make_region=True,region_file=None,
                      galactic_index_free=False,write_directory='',use_old_names=False):
                
		self.free_radius=(self.roi[2] if free_radius<=0 else free_radius)
		self.max_free_radius=(self.free_radius if max_free_radius is None else max_free_radius)
		if self.max_free_radius<self.free_radius:
                        warnings.warn(f'max_free_radius ({max_free_radius:.1f}) is < free_radius ({free_radius:.1f}), making the parameter useless.')
                self.extra_radius=extra_radius

		self.variable_free=variable_free
		self.force_point_sources=force_point_sources
		self.extended_catalog_names=extended_catalog_names
		self.norms_free_only=norms_free_only

		self.extended_directory=extended_directory if extended_directory is None else\
                                    os.sep.join(['$(FERMI_DIR)','data','pyBurstAnalysisGUI','templates'])
		
		self.sigma_to_free=sigma_to_free
		self.make_region=make_region
		self.galactic_index_free=galactic_index_free

                self.create_region_model()
		
		if make_region:
                        if region_file is None:
                                rhold=os.path.basename(self.output_name).split('.')[:-1]
                                rhold=reduce(lambda s1,s2:s1+s2,rhold).join(['ROI_','.reg'])
                                
                                self.region_file=os.path.join(write_dir,rhold)
                        else:
                                self.region_file=region_file
                        
                        print('Building ds9-style region file',end='...')
                        build_region(self.region_file,self.final_sources)
                        print(f'done!\nFile saved as {self.region_file}.')
        
        def create_region_model(self):
		print(f'Creating spatial and spectral model from the 4FGL DR-{self.DR} catalog: {self.catalog_file}.')
		
		if os.path.basename(self.catalog_file).split('.')[-1].lower()=='xml':
                        self.get_sources_xml()
		else:
                        self.get_sources_fits()

                #from either method we get a self.sources attribute which is a nested dictionary with all the sources
                self.build_model()

        def build_model(self):
                #use the self.sources nested dictionary to create the output xml file
	

        def get_sources_xml(self):
                catalog=minidom.parse(self.catalog_file)
                sources=catalog.getElementsByTagName('source')

                #traverse through the sources once, just to get positional information
                right_ascensions,declinations=[],[]
                for source in sources:
                        try:
                                right_ascentions.append(float(source.getAttribute('RA')))
                                declinations.append(float(source.getAttribute('DEC')))
                        except:
                                for parameter in source.getElementsByTagName('spatialModel')[0].getElementsByTagName('parameter'):
                                        if parameter.getAttribute('name')=='RA':
                                                right_ascensions.append(float(parameter.getAttribute('value')))
                                        elif parameter.getAttribute('name')=='DEC':
                                                declinations.append(float(parameter.getAttribute('value')))

                right_ascensions=np.array(right_ascensions)
                declinations=np.array(declinations)
                
                self.source_distances=angular_separation(self.ROI[0],self.ROI[1],
                                                         right_ascensions,declinations)

                source_mask=self.source_distances<=self.ROI[2]+self.extra_radius
                self.source_distances=self.source_distances[source_mask]
                distance_index=self.source_distances.argsort()
                self.source_distances=self.source_distances[distance_index]
                
                sources=sources[source_mask][distance_index]
                right_ascensions=right_ascentions[source_mask][distance_index]
                declinations=declinations[source_mask][distance_index]

                #now traverse the shorter list of sources to create the self.sources nested dictionary
                self.sources={}
                for source,RA,DEC,distance in zip(self.sources,right_ascensions,declinations,self.source_distances):
                        name=source.getAttribute('name')
                        self.sources.update([(name,{'roi_distance':distance})])

                        if source.getAttribute('type')=='DiffuseSource':
                                self.sources[name].update([('Extended',not self.force_point_sources)])

                        #now get the spectrum information
                        spectrum=source.getElementsByTagName('spectrum')
                        
                        self.sources[name].update([('spectrum',
                             {'model':specrum[0].getAttribute('type')})])
                        
                        for parameter in spectrum[0].getElementsByTagName('parameter'):
                        	self.sources[name]['spectrum'].update([(parameter.getAttribute('name'),
                        	parameter.getAttribute('value')*parameter.getAttribute('scale'))])
                        
                        #now for the spatial information
                        spatial=source.getElementsByTagName('spatialModel')
                        
                        if self.sources[name]['Extended']:
                        	self.sources[name].update([('spatial',
                                        {'spatial_model':spatial[0].getAttribute('type')})])
                        	for parameter in spatial[0].getElementsByTagName('parameter'):
                        		self.sources[name]['spatial'].update([(parameter.getAttribute('name'),
                        		parameter.getAttribute('value'))])
                        
                        else:
                        	self.sources[name].update([('spatial',
                        	    {'spatial_model':'SkyDir','RA':RA,'DEC':DEC})])
                        	

                        #now just need to add the TS and variability index stuff (if possible)
                        

        def get_sources_fits(self):
                #first, open the catalog FITS file and extract the information
                #we care about, cast as numpy arrays for ease later
                with pyfits.open(self.catalog_file) as catalog:
                        source_info=catalog['LAT_Point_Source_Catalog'].data.field
                        extended_info=catalog['ExtendedSources'].data.field
                        
                        extended_names=np.array(extended_info('Source_Name'))
                        extended_files=np.array(extended_info('Spatial_Filename'))
                	extended_spatial_functions=np.array(extended_info('Spatial_Function'))
                	extended_extents=np.array(extended_info('Model_SemiMajor'))
                        extended_RA=np.array(extended_info('RAJ2000'))
                	extended_DEC=np.array(extended_info('DEJ2000'))

                        #cast these as a data frame with the extended_name as the index
                	#making it much easier to search later using just .loc
                	#instead of a nested for loop
                	extended_sources=pd.DataFrame(np.c_[extended_files,extended_spatial_functions,
                                extended_extents,extended_RA,extended_DEC],
                                columns=['file','function','extent','RA','DEC'],
                                index=extended_names)
                	
                	names=np.array(source_info('Source_Name'))
                	
                	average_significances=np.array(source_info('Signif_Avg'))
                	variability_indices=np.array(source_info('Variability_Index'))
                	
                	extended_source_names=np.array(source_info('Extended_Source_Name'))
                	
                	right_ascensions=np.array(source_info('RAJ2000'))
                	declinations=np.array(source_info('DEJ2000'))
                	
                	power_law_fluxes=np.array(source_info('PL_Flux_Density'))
                	log_parabola_fluxes=np.array(source_info('LP_Flux_Density'))
                	cutoff_fluxes=np.array(source_info('PLEC_Flux_Density'))
                	
                	pivot_energies=np.array(source_info('Pivot_Energy'))
                	
                	power_law_indices=np.array(source_info('PL_Index'))
                	
                	log_parabola_indices=np.array(source_info('LP_Index'))
                	log_parabola_betas=np.array(source_info('LP_beta'))
                	
                        cutoff_indices=np.array(source_info(('PLEC_IndexS' if self.DR==3 else 'PLEC_Index')))
                        cutoff_expfactors=np.array(source_info(('PLEC_ExpfactorS' if self.DR==3 else 'PLEC_Expfactor')))
                        cutoff_expindices=np.array(source_info('PLEC_Exp_Index'))
                        
                        spectral_types=np.array(source_info('SpectrumType'))

                #calculate the source distances with respect to the ROI center
                self.source_distances=angular_separation(self.ROI[0],self.ROI[1],
                                                         right_ascensions,declinations)

                #get a mask for sources in the distance range we care about
                #and then sort them in ascending order
                source_mask=self.source_distances<=self.ROI[2]+self.extra_radius
                self.source_distances=self.source_distances[source_mask]
                distance_index=self.source_distances.argsort()
                self.source_distances=self.source_distances[distance_index]

                #now, downselect all the other arrays
                right_ascensions=right_ascentions[source_mask][distance_index]
                declinations=declinations[source_mask][distance_index]
                
                names=names[source_mask][distance_index]
                
                average_significances=average_significances[source_mask][distance_index]
                variability_indices=variability_indices[source_mask][distance_index]
                
                extended_source_names=extended_source_names[source_mask][distance_index]
                
                power_law_fluxes=power_law_fluxes[source_mask][distance_index]
                log_parabola_fluxes=log_parabola_fluxes[source_mask][distance_index]
                cutoff_fluxes=cutoff_fluxes[source_mask][distance_index]
                
                pivot_energies=pivot_energies[source_mask][distance_index]
                
                power_law_indices=power_law_indices[source_mask][distance_index]
                
                log_parabola_indices=log_parabola_indices[source_mask][distance_index]
                log_parabola_betas=log_parabola_betas[source_mask][distance_index]
                
                cutoff_indices=cutoff_indices[source_mask][distance_index]
                cutoff_expfactors=cutoff_expfactors[source_mask][distance_index]
                cutoff_expindices=cutoff_expindices[source_mask][distance_index]
                
                spectral_types=spectral_types[source_mask][distance_index]

                self.sources={}

                for idx in range(len(names)):
                        #check if the source is extended or not, decide on name for source
                        #and if it is extended
                        if names[idx][-1]=='e':
                                name=(names[idx] if self.extended_catalog_names else extended_source_names[idx])
                                self.sources.update([(name,{'roi_distance':self.source_distances[idx]})])
                                self.sources[name].update([('Extended',not self.force_point_sources)])
                        else:
                                name=names[idx]
                                self.sources.update([(name,{'roi_distance':self.source_distances[idx]})])
                                self.sources[name].update([('Extended',False)])

                        #now get source spectral information in the form we need
                        #I need to verify that PowerLaw, LogParabola, and PLSuperExpCutoff4
                        #are the only models we have in 4FGL DR 1, 2, and 3
                        if spectral_types[idx]=='PowerLaw':
                                self.sources[name].update([('spectrum',
                                        {'model':'PowerLaw',
                                         'Prefactor':power_law_fluxes[idx],
                                         'Index':power_law_indices[idx],
                                         'Scale':pivot_energies[idx],
                                         'variable':variability_indices[idx]>=self.variability_threshold,
                                         'significant':average_significances[idx]>=self.sigma_to_free})])

                        elif spectral_types[idx]=='LogParabola':
                                self.sources[name].update([('spectrum',
                                        {'model':'LogParabola',
                                         'norm':log_parabola_fluxes[idx],
                                         'alpha':log_prabola_indices[idx],
                                         'beta':log_parabola_betas[idx],
                                         'Eb':pivot_energies[idx],
                                         'variable':variability_indices[idx]>=self.variability_threshold,
                                         'significant':average_significances[idx]>=self.sigma_to_free})])

                        else:
                                if self.DR==3:
                                        self.sources[name].update([('spectrum',
                                                {'model':'PLSuperExpCutoff4',
                                                 'Prefactor':cutoff_fluxes[idx],
                                                 'IndexS':cutoff_indices[idx],
                                                 'Scale':pivot_energies[idx],
                                                 'ExpfactorS':cutoff_expfactors[idx],
                                                 'Index2':cutoff_expindices[idx],
                                                 'variable':variability_indices[idx]>=self.variability_threshold,
                                                 'significant':average_significances[idx]>=self.sigma_to_free})])

                        #now for the spatial information
                        #check if extended, this flag already uses force_point_sources flag
                        if self.sources[name]['Extended']:
                                #do extended stuff
                                row=extended_sources.loc[extended_source_names[idx]]
                                #if one of the 'Radial' models
                                if row.function[:6]=='Radial':
                                        self.sources[name].update([('spatial',
                                                {'spatial_model':'RadialGaussian' if row.function=='RadialGauss' else row.function,
                                                 'RA':row.RA,
                                                 'DEC':row.DEC})])
                                        if row.function=='RadialDisk':
                                                self.sources[name]['spatial'].update([('Radius',row.extent)])
                                        else:
                                                self.sources[name]['spatial'].update([('Sigma',row.extent/(-2.np.log(0.32))**0.5)])

                                #otherwise it is a spatial map
                                else:
                                        self.sources[name].update([('spatial',
                                                {'type':'SpatialMap',
                                                 'spatial_file':row.file})])

                        #point source
                        else:
                                self.sources[name].update([('spatial',
                                        {'type':'SkyDir',
                                         'RA':right_ascensions[idx],
                                         'DEC':declinations[idx]})])

                #done with the loop, object will use the sources dictionary to build the model



                        

def addSrcsXML(sL,GD,GDn,ISO,ISOn,use_old_names=False):
	varValue=sL.variability_threshold
	inputXml=minidom.parse(sL.srcs)
	outputXml=minidom.getDOMImplementation().createDocument(None,'source_library',None)
	outputXml.documentElement.setAttribute('title','source library')
	catalog=inputXml.getElementsByTagName('source')
	Sources={}
	ptSrcNum=0
	extSrcNum=0
	for src in catalog:
		if src.getAttribute('type')=='PointSource':
			for p in src.getElementsByTagName('spatialModel')[0].getElementsByTagName('parameter'):
				if p.getAttribute('name')=='RA':
					srcRA=float(p.getAttribute('value'))
				if p.getAttribute('name')=='DEC':
					srcDEC=float(p.getAttribute('value'))
		else:
			try:
			  srcDEC=float(src.getAttribute('DEC'))
			  srcRA=float(src.getAttribute('RA'))
			except:
			  for p in src.getElementsByTagName('spatialModel')[0].getElementsByTagName('parameter'):
			    if p.getAttribute('name')=='RA':
			      srcRA=float(p.getAttribute('value'))
			    if p.getAttribute('name')=='DEC':
			      srcDEC=float(p.getAttribute('value'))
		dist=angsep(sL.roi[0],sL.roi[1],srcRA,srcDEC) #check that source is within ROI radius + 10 degress of ROI center
		if srcRA==sL.roi[0] and srcDEC==sL.roi[1]:
			dist=0.0
		if dist<=sL.roi[2]+sL.ER:
			spec=src.getElementsByTagName('spectrum')
			specType=spec[0].getAttribute('type')
			specPars=spec[0].getElementsByTagName('parameter')
			Ext=(True if (src.getAttribute('type')=='DiffuseSource' and not sL.psF) else False)
			sname=src.getAttribute('name')
			
			if use_old_names:#if you want the same naming convention as in make1FGLxml.py and make2FGLxml.py, e.g., preceeded by an underscore and no spaces
				sn='_'
				for N in str(sname).split(' '):
					sn+=N
			varIdx=float(src.getAttribute('Variability_Index'))
			#varidx=-1#hold over from old XMLs having and then not having and then now having Variability Index
			Sources[sname]={'ra':srcRA,'dec':srcDEC,'E':Ext,'stype':str(specType)}
			specOut=outputXml.createElement('spectrum')
			if str(specType)=='PLSuperExpCutoff2':
			  specOut.setAttribute('type','PLSuperExpCutoff')
			else:
			  specOut.setAttribute('type',specType)
			spatialOut=outputXml.createElement('spatialModel')
			srcOut=outputXml.createElement('source')
			srcOut.setAttribute('name',sname)
			srcOut.setAttribute('ROI_Center_Distance',"%.2f"%dist)
			if dist>=sL.roi[2] or dist>=sL.max_free_radius:
				Sources[sname]['free']=False
				for p in specPars:
				  specOut.appendChild(parameter_element("0","%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
			elif dist>sL.free_radius:
				if sL.var and varIdx>=varValue:
					Sources[sname]['free']=True
					for p in specPars:
					  freeFlag=("1" if p.getAttribute('name')==spec[0].getAttribute('normPar') else "0")
					  specOut.appendChild(parameter_element("%s"%freeFlag,"%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
				else:
					Sources[sname]['free']=False
					for p in specPars:
					  specOut.appendChild(parameter_element("0","%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
			elif float(src.getAttribute('TS_value'))>=sL.sig:
				Sources[sname]['free']=True
				#specOut.setAttribute('apply_edisp',ed)
				for p in specPars:
				  freeFlag=("1" if p.getAttribute('name')==spec[0].getAttribute('normPar') or (not sL.nO and p.getAttribute('free')=="1") else "0")
				  specOut.appendChild(parameter_element("%s"%freeFlag,"%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
			else:
				if sL.var and varIdx>=varValue:
					Sources[sname]['free']=True
					#specOut.setAttribute('apply_edisp',ed)
					for p in specPars:
					  freeFlag=("1" if p.getAttribute('name')==spec[0].getAttribute('normPar') else "0")
					  specOut.appendChild(parameter_element("%s"%freeFlag,"%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
				else:
					Sources[sname]['free']=False
					#specOut.setAttribute('apply_edisp','false')
					for p in specPars:
					  specOut.appendChild(parameter_element("0","%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
			if Ext:
				spatial=src.getElementsByTagName('spatialModel')
				spatType=spatial[0].getAttribute('type')
				spatPars=spatial[0].getElementsByTagName('parameter')
				if str(spatType)=='SpatialMap':
				  spatialOut.setAttribute('type','SpatialMap')
				  spatialOut.setAttribute('map_based_integral','true')
				  efile=sL.extD+spatial[0].getAttribute('file').split('/')[-1]
				  spatialOut.setAttribute('file',efile)
				  print('Extended source %s in ROI, make sure %s is the correct path to the extended template.'%(sname,efile))
				else:#have to do above to get correct extended source template file localtion
				  spatialOut.setAttribute('type',str(spatType))
				  for p in spatPars:#for radial disks and gaussians, can just do the following
				    spatialOut.appendChild(parameter_element("0","%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
				    print('Extended source %s in ROI, with %s spatial model.'%(sname,str(spatType)))
				srcOut.setAttribute('type','DiffuseSource')
				extSrcNum+=1
				#print('Extended source %s in ROI, make sure %s is the correct path to the extended template.'%(sname,efile))
			else:
				spatialOut.setAttribute('type','SkyDirFunction')
				spatialOut.appendChild(parameter_element("0","RA","360.0","-360.0","1.0","%.4f"%srcRA))
				spatialOut.appendChild(parameter_element("0","DEC","360.0","-360.0","1.0","%.4f"%srcDEC))
				srcOut.setAttribute('type','PointSource')
				ptSrcNum+=1
			srcOut.appendChild(specOut)
			srcOut.appendChild(spatialOut)
			outputXml.documentElement.appendChild(srcOut)
	gal=outputXml.createElement('source')
	gal.setAttribute('name',GDn)
	gal.setAttribute('type','DiffuseSource')
	galspec=outputXml.createElement('spectrum')
	galspec.setAttribute('type','PowerLaw')
	#galspec.setAttribute('apply_edisp','false')
	galspec.appendChild(parameter_element("1","Prefactor","10","0","1","1"))
	if sL.GIF:
		galspec.appendChild(parameter_element("1","Index","1","-1","1","0"))
	else:
		galspec.appendChild(parameter_element("0","Index","1","-1","1","0"))
	galspec.appendChild(parameter_element("0","Scale","1e6","2e1","1","100"))
	galspatial=outputXml.createElement('spatialModel')
	galspatial.setAttribute('type','MapCubeFunction')
	galspatial.setAttribute('file',GD)
	galspatial.appendChild(parameter_element("0","Normalization","1e3","1e-3","1","1"))
	gal.appendChild(galspec)
	gal.appendChild(galspatial)
	outputXml.documentElement.appendChild(gal)
	iso=outputXml.createElement('source')
	iso.setAttribute('name',ISOn)
	iso.setAttribute('type','DiffuseSource')
	isospec=outputXml.createElement('spectrum')
	isospec.setAttribute('type','FileFunction')
	isospec.setAttribute('file',ISO)
	isospec.setAttribute('apply_edisp','false')
	isospec.appendChild(parameter_element("1","Normalization","10","0.01","1","1"))
	isospatial=outputXml.createElement('spatialModel')
	isospatial.setAttribute('type','ConstantValue')
	isospatial.appendChild(parameter_element("0","Value","10","0","1","1"))
	iso.appendChild(isospec)
	iso.appendChild(isospatial)
	outputXml.documentElement.appendChild(iso)
	xmlStr=outputXml.toprettyxml(' ').splitlines(True)
	outStr=filter(lambda xmlStr: len(xmlStr) and not xmlStr.isspace(),xmlStr)
	outfile=open(sL.out,'w')
	outfile.write(''.join(outStr))
	outfile.close()
	if not sL.psF:
		print('Added %i point sources and %i extended sources'%(ptSrcNum,extSrcNum))
		if extSrcNum>0:
			print('If using unbinned likelihood you will need to rerun gtdiffrsp for the extended sources or rerun the makeModel function with optional argument force_point_sources=True')
	else:
		print('Added %i point sources, note that any extended sources in ROI were modeled as point sources becaue force_point_sources option was set to True'%ptSrcNum)
	if sL.reg:
		BuildRegion(sL,Sources)
	return

#function to cycle through the source list and add point source entries
def addSrcsFITS(sL,GD,GDn,ISO,ISOn,use_old_names):
	model=open(sL.out,'w') #open file in write mode, overwrites other files of same name
	catfile=pyfits.open(sL.srcs) #open source list file and access necessary fields, requires LAT source catalog definitions and names
	data=catfile['LAT_Point_Source_Catalog'].data
	extendedinfo=catfile['ExtendedSources'].data
	extName=extendedinfo.field('Source_Name')	
	extFile=extendedinfo.field('Spatial_Filename')
	extFunc=extendedinfo.field('Spatial_Function')
	extSize=extendedinfo.field('Model_SemiMajor')
	extRa=extendedinfo.field('RAJ2000')
	extDec=extendedinfo.field('DEJ2000')
	name=data.field('Source_Name')
	Sigvals=data.field('Signif_Avg')
	VarIdx=data.field('Variability_Index')
	EName=data.field('Extended_Source_Name')
	ra=data.field('RAJ2000')
	dec=data.field('DEJ2000')
	plflux=data.field('PL_Flux_Density')
	lpflux=data.field('LP_Flux_Density')
	coflux=data.field('PLEC_Flux_Density')
	pivot=data.field('Pivot_Energy')
	plIndex=data.field('PL_Index')
	lpIndex=data.field('LP_Index')
	lpbeta=data.field('LP_beta')
	if sL.DR==3:
		plecIndex=data.field('PLEC_IndexS')
		plecexpFact=data.field('PLEC_ExpfactorS')
		plecexpIndex=data.field('PLEC_Exp_Index')
	else:
		plecIndex=data.field('PLEC_Index')
		plecexpFact=data.field('PLEC_Expfactor')
		plecexpIndex=data.field('PLEC_Exp_Index')
	spectype=data.field('SpectrumType')
	model.write('<?xml version="1.0" ?>\n')
	model.write('<source_library title="source library">\n')
	model.write('\n<!-- Point Sources -->\n')
	step=(sL.roi[2]+sL.ER)/5. #divide ROI radius plus extra_radiusius degrees into 5 steps for ordering of sources
	i=1
	radii=[]
	ptSrcNum=0
	extSrcNum=0
	Sources={}#dictionary for sources, useful for creating region file later.
	while i<6:
		if i*step<=sL.roi[2]+sL.ER:
			radii+=[step*i]
		else:
			radii+=[sL.roi[2]+sL.ER] #just in case of rounding errors
		i+=1
	for x in radii:
		if x==sL.roi[2]+sL.ER:
			model.write('\n<!-- Sources between [%s,%s] degrees of ROI center -->\n' %(x-step,x))
		else:
			model.write('\n<!-- Sources between [%s,%s) degrees of ROI center -->\n' %(x-step,x))
		for n,plf,lpf,cof,r,d,p,pli,lpi,lpb,pleci,plecef,plecei,t,TS,En,vi in zip(name,plflux,lpflux,coflux,ra,dec,pivot,plIndex,lpIndex,lpbeta,plecIndex,plecexpFact,plecexpIndex,spectype,Sigvals,EName,VarIdx):
			E=(True if n[-1]=='e' else False)
			dist=angsep(sL.roi[0],sL.roi[1],r,d) #check that source is within ROI radius + 10 degress of ROI center
			if r==sL.roi[0] and d==sL.roi[1]:
				dist=0.0
			if (dist<x and dist>=x-step) or (x==sL.roi[2]+10. and dist==x):
				if E and not sL.psF:
					Sources[En]={'ra':r,'dec':d,'stype':t,'E':E}
					extSrcNum+=1
					Name='<source ROI_Center_Distance="%.3f" name="%s" type="DiffuseSource">\n' %(dist,En)
				else:
					if E and not sL.E2C:#even if forcing all to point sources, use extended name except if extended_catalog_names flag is set
						Sources[En]={'ra':r,'dec':d,'stype':t,'E':E}
						Name='<source ROI_Center_Distance="%.3f" name="%s" type="PointSource">\n' %(dist,En)
					else:
						Sources[n]={'ra':r,'dec':d,'stype':t,'E':E}
						if use_old_names:
							srcname='_'
							for N in n.split(' '):
								srcname+=N
							Name='<source ROI_Center_Distance="%.3f" name="%s" type="PointSource">\n' %(dist,srcname)
						else:
							Name='<source ROI_Center_Distance="%.3f" name="%s" type="PointSource">\n' %(dist,n)
					ptSrcNum+=1
				if t=='PowerLaw':
					#uncomment out the two lines immediately following later
					spec,free=PLspec(sL,plf,pli,p,dist,TS,vi)
				elif t=='LogParabola':
					fixAll=(True if n=='4FGL J0534.5+2201i' else False)
					spec,free=LPspec(sL,lpf,lpi,p,lpb,dist,TS,vi,fixAll)
				else:
					if sL.DR==3:
						spec,free=CO4spec(sL,cof,pleci,p,plecef,plecei,dist,TS,vi)
					else:
						spec,free=CO2spec(sL,cof,pleci,p,plecef,plecei,dist,TS,vi)
				if E and not sL.E2C:
					Sources[En]['free']=free
				else:
					Sources[n]['free']=free
				if E and not sL.psF:
					efile=None
					efunc=None
					eSize=None
					eR=None
					eD=None
					for EXTNAME,EXTFILE,EXTFUNC,EXTSIZE,EXTRA,EXTDEC in zip(extName,extFile,extFunc,extSize,extRa,extDec):
						if En==EXTNAME:
							efunc=EXTFUNC
							efunc=('RadialGaussian' if efunc=='RadialGauss' else efunc)
							if efunc=='SpatialMap':
							  efile=sL.extD+EXTFILE
							else:
							  eSize=(EXTSIZE/(-2*log(0.32))**0.5 if efunc=='RadialGaussian' else EXTSIZE)
							  eR=EXTRA
							  eD=EXTDEC
					if efunc=='SpatialMap':
					  if efile==None:
					    print('could not find a match for',En,'in the list:')
					    print(extName)
					    efile=''
					  skydir='\t<spatialModel file="%s" map_based_integral="true" type="SpatialMap">\n'%(efile)
					  print('Extended source %s in ROI, make sure %s is the correct path to the extended template.'%(En,efile))
					  skydir+='\t\t<parameter free="0" max="1000" min="0.001" name="Prefactor" scale="1" value="1"/>\n'
					  skydir+='\t</spatialModel>\n'
					else:
					  skydir='\t<spatialModel type="%s">\n'%efunc
					  skydir+='\t<parameter free="0" max="360" min="-360" name="RA" scale="1" value="%s"/>\n'%eR
					  skydir+='\t<parameter free="0" max="90" min="-90" name="DEC" scale="1" value="%s"/>\n'%eD
					  if efunc=='RadialDisk':
					    skydir+='\t<parameter free="0" max="10" min="0" name="Radius" scale="1" value="%s"/>\n'%eSize
					  else:
					    skydir+='\t<parameter free="0" max="10" min="0" name="Sigma" scale="1" value="%s"/>\n'%eSize
					  skydir+='\t</spatialModel>\n'
					  print('Extended source %s in ROI with %s spatial model.'%(En,efunc))
				else:
					skydir='\t<spatialModel type="SkyDirFunction">\n'
					skydir+='\t\t<parameter free="0" max="360.0" min="-360.0" name="RA" scale="1.0" value="%s"/>\n' %r
					skydir+='\t\t<parameter free="0" max="90" min="-90" name="DEC" scale="1.0" value="%s"/>\n' %d
					skydir+='\t</spatialModel>\n'
				skydir+='</source>'
				(src,)=(Name+spec+skydir,)
				ptsrc=pS(src).getElementsByTagName('source')[0]
				ptsrc.writexml(model)
				model.write('\n')
	catfile.close() #close file
	if not sL.psF:
		print('Added %i point sources and %i extended sources'%(ptSrcNum,extSrcNum))
		if extSrcNum>0:
			print('If using unbinned likelihood you will need to rerun gtdiffrsp for the extended sources or rerun the makeModel function with optional argument force_point_sources=True')
	else:
		print('Added %i point sources, note that any extended sources in ROI were modeled as point sources becaue force_point_sources option was set to True'%ptSrcNum)
	#add galactic diffuse with PL spectrum, fix index to zero for general use, those who want it to be free can unfreeze parameter manually
	model.write('\n<!-- Diffuse Sources -->\n')
	Name='\n<source name="%s" type="DiffuseSource">\n' %GDn
	spec='\t<spectrum type="PowerLaw">\n'
	spec+='\t\t<parameter free="1" max="10" min="0" name="Prefactor" scale="1" value="1"/>\n'
	if sL.GIF:
		spec+='\t\t<parameter free="1" max="1" min="-1" name="Index" scale="1.0" value="0"/>\n'
	else:
		spec+='\t\t<parameter free="0" max="1" min="-1" name="Index" scale="1.0" value="0"/>\n'
	spec+='\t\t<parameter free="0" max="2e2" min="5e1" name="Scale" scale="1.0" value="1e2"/>\n'
	spec+='\t</spectrum>\n'
	skydir='\t<spatialModel file="%s" type="MapCubeFunction">\n' %GD
	skydir+='\t\t<parameter free="0" max="1e3" min="1e-3" name="Normalization" scale="1.0" value="1.0"/>\n'
	skydir+='\t</spatialModel>\n'
	skydir+='</source>'
	(src,)=(Name+spec+skydir,)
	galdiff=pS(src).getElementsByTagName('source')[0]
	galdiff.writexml(model)
	model.write('\n')
	Name='<source name="%s" type="DiffuseSource">\n' %ISOn
	spec='\t<spectrum type="FileFunction" file="%s"  apply_edisp="false">\n' %ISO
	spec+='\t\t<parameter free="1" max="10" min="1e-2" name="Normalization" scale="1" value="1"/>\n'
	spec+='\t</spectrum>\n'
	skydir='\t<spatialModel type="ConstantValue">\n'
	skydir+='\t\t<parameter free="0" max="10.0" min="0.0" name="Value" scale="1.0" value="1.0"/>\n'
	skydir+='\t</spatialModel>\n'
	skydir+='</source>'
	(src,)=(Name+spec+skydir,)
	iso=pS(src).getElementsByTagName('source')[0]
	iso.writexml(model)	
	model.write('\n</source_library>')
	model.close()
	if sL.reg:
		BuildRegion(sL,Sources)
	return

def mybool(Input):
	return {'True':True,'False':False,'T':True,'F':False,'t':True,'f':False,'TRUE':True,'FALSE':False,"true":True,"false":False,"1":True,"0":False}.get(Input)

def cli():
	import argparse
	
	helpString="Creates an xml model from the 3FGL catalog (FITS or xml version) for a specific ROI,\
		    coordinates of the ROI center are taken from an input event file,\
		    the radius for including sources is 10 degrees beyond the extraction radius used in the event file,\
		    sources with free parameters within the original extraction radius are chosen based on nearness to center, significance, and variability."
	parser=argparse.ArgumentParser(description=helpString)
	parser.add_argument("catalog",type=str,help="Catalog file to use, can be FITS or xml.")
	parser.add_argument("ev",type=str,help="Event file with ROI information in header.")
	parser.add_argument("-o","--outputxml",type=str,default='mymodel.xml',help="Name of output xml file, is set to overwrite files of same name.")
	parser.add_argument("-G","--galfile",type=str,default='$(FERMI_DIR)/refdata/fermi/galdiffuse/gll_iem_v07.fits',help="Name and location of Galactic diffuse model to use, will default to gll_iem_v06.fits.")
	parser.add_argument("-g","--galname",type=str,default='gll_iem_v07',help="Name of Galactic diffuse component in output model, will default to gll_iem_v06.")
	parser.add_argument("-I","--isofile",type=str,default='$(FERMI_DIR)/refdata/fermi/galdiffuse/iso_P8R3_SOURCE_V3_v1.txt',help="Name of isotropic diffuse template for output model, will default to P8R3 SOURCE class model.")
	parser.add_argument("-i","--isoname",type=str,default='iso_P8R3_SOURCE_V3_v1',help="Name of isotropic diffuse component in output model, default is for P8R3 SOURCE class.")
	parser.add_argument("-N","--normsonly",type=mybool,default=False,help="Flag to only let the normalizations of parameters be free, default is False.",nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument("-e","--extended_directory",type=str,default='',help="Path to directory with LAT extended source templates, will default to STs default.")#need to figure out what that is
	parser.add_argument("-r","--free_radius",type=float,default=-1.,help="Radius, in degrees, from ROI center beyond which all source parameters should be fixed, will default to selection radius.")
	parser.add_argument("-R","--max_free_radius",type=float,default=None,help="Absolute maximum radius, in degrees, from ROI center beyond which all source parameters should be fixed, even variable sources will not be freed beyond this radius, defaults to free_radius value.")
	parser.add_argument("-ER","--extra_radius",type=float,default=10.,help="Radius beyond event file ROI out to which sources will be included in the model with all parameters fixed, default is 10, good for analyses starting around a few hundred MeV, can be decreased for high energy only fits.")
	parser.add_argument("-s","--sigma_to_free",type=float,default=5.,help="Average significance below which all source parameters are fixed, defaults to 5.  Note, if using the 3FGL catalog xml file as input, this is actually a cut on TS, so adjust accordingly.")
	parser.add_argument("-v","--variable_free",type=mybool,default=True,help="Flag to set normalization of significantly variable sources, even if source is beyond radius limit or below TS limit, default is True.",choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument("-p","--force_point_sources",type=mybool,default=False,help="Flag to cast extended sources as point sources, default is False.",nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument("-E2C","--extended_catalog_names",type=mybool,default=False,help="Flag to use catalog names for extended sources, default is False.",nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument("-m","--make_region",type=mybool,default=True,help="Flag to create ds9 region file as well as the xml model, default is True.",choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument("-GIF","--galactic_index_free",type=mybool,default=False,help="Flag to use a power-law modification to the Galactic diffuse model spectrum and have the index be free, default is False.",nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument("-wd","--writeDir",type=str,default='',help="Directory to write the output ds9 region file in if not the current working directory or if you are specifying the full path to the newly made XML file.")
	parser.add_argument("-ON","--use_old_names",type=mybool,default=False,help="Flag to use the make2FLGxml style naming convention, underscore before name and no spaces, default is False.",nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
	parser.add_argument('-DR','--DataRelease',type=int,default=3,help='Choice of data release, 3 for 4FGL-DR3, 2 for 4FGL-DR2, and 1 for 4FGL.',choices=[1,2,3])
	
	args=parser.parse_args()

	sL=srcList(args.catalog,args.ev,args.outputxml,args.DataRelease)
	
	sL.makeModel(args.galfile,args.galname,args.isofile,args.isoname,args.normsonly,args.extended_directory,args.free_radius,args.max_free_radius,args.extra_radius,args.sigma_to_free,args.variable_free,args.force_point_sources,args.extended_catalog_names,args.make_region,args.galactic_index_free,args.writeDir,args.use_old_names)
	
	
if __name__=='__main__': cli()
