from xml.dom import minidom
from xml.minidom import parseString
import numpy as np

#copied from Damien's macro
#Create an XML document parameter description element
def parameter_element(free,name,maximum,minimum,scale,value):
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

###########################################################################################
#create separate classes for the different spectral models
#all built on a base class to avoid too much duplication of functions/code
#only need to pass the spectrum attribute to the appendChild method of a source XML element
###########################################################################################
class Spectrum:
    def __init__(self):
        self.document=minidom.getDOMImplementation().createDocument(None,None,None)
        self.spectrum=self.document.createElement('spectrum')

    def build(self,model_name,parameters):
        self.spectrum.setAttribute('type',model_name)
        for parameter in parameters:
            self.spectrum.appendChild(parameter_element(parameter.get('free'),
                parameter.get('name'),parameter.get('max'),parameter.get('min'),
                parameter.get('scale'),parameter.get('value')))

#PowerLaw class
class PowerLaw(Spectrum):
    def __init__(self,prefactor=1e-12,index=2,scale_energy=1000):
        Spectrum.__init__()

        #some checks on the inputs
        if prefactor<0 or scale_energy<=0:
            raise ValueError("Input parameters invalid, at least one of prefactor or scale_energy is negative or zero.")
        
        self.prefactor_scale=10**np.floor(np.log10(prefactor))
        self.prefactor=prefactor/self.prefactor_scale
        self.index=abs(index)
        self.scale_energy=scale_energy

        self.model_name='PowerLaw'
        self.norm_par='Prefactor'

    def build_spectrum(self,prefactor_free=True,index_free=True):
        self.parameters=[{'free':int(prefactor_free),'name':'Prefactor','min':0,
                     'max':1e6,'scale':self.prefactor_scale,'value':self.prefactor}]

        self.parameters.append({'free':int(index_free),'name':'Index','min':0,
                     'max':10,'scale':-1,'value':self.index})

        self.parameters.append({'free':0,'name':'Scale','min':min(100,self.scale_energy*0.9),
                           'max':max(500000,self.scale_energy*1.5),'scale':1,'value':self.scale_energy})

        self.build(self.model_name,self.parameters)

#PowerLaw2 class
class PowerLaw2(Spectrum):
    def __init__(self,integral=1e-8,index=2,lower_limit=100,upper_limit=300000):
        Spectrum.__init__()

        #some checks on the inputs
        if integral<0 or lower_limit<=0 or upper_limit<=0:
            raise ValueError("Input parameters invalid, at least one of integral, lower_limit, or upper_limit is negative or zero.")
        
        if lower_limit>=upper_limit:
            raise ValueError(f"Requested LowerLimit of {lower_limit:,} is >= UpperLimit of {upper_limit:,}.")
        
        self.integral_scale=10**np.floor(np.log10(integral))
        self.integral=prefactor/self.integral_scale
        self.index=abs(index)
        self.lower_limit=lower_limit
        self.upper_limit=upper_limit

        self.model_name='PowerLaw2'
        self.norm_par='Integral'

    def build_spectrum(self,integral_free=True,index_free=True):
        self.parameters=[{'free':int(integral_free),'name':'Integral','min':0,
                     'max':1e6,'scale':self.prefactor_scale,'value':self.integral}]

        self.parameters.append({'free':int(index_free),'name':'Index','min':0,
                     'max':10,'scale':-1,'value':self.index})

        self.parameters.append({'free':0,'name':'LowerLimit','min':10,
                           'max':max(500000,self.lower_limit*1.5),
                           'scale':1,'value':self.lower_limit})

        self.parameters.append({'free':0,'name':'UpperLimit','min':10,
                           'max':max(500000,self.upper_limit*1.5),
                           'scale':1,'value':self.upper_limit})

        self.build(self.model_name,self.parameters)

#PLSuperExpCutoff class
class PLSuperExpCutoff(Spectrum):
    def __init__(self,prefactor=1e-12,index1=2,scale_energy=1000,cutoff=1000,index2=1):
        Spectrum.__init__()

        #some checks on the inputs
        if prefactor<0 or scale_energy<=0 or cutoff<=0 or index2<0:
            raise ValueError("Input parameters invalid, at least one of prefactor, scale_energy, cutoff, or index2 is negative or zero.")
        
        self.prefactor_scale=10**np.floor(np.log10(prefactor))
        self.prefactor=prefactor/self.prefactor_scale
        self.index1=abs(index1)
        self.scale_energy=scale_energy
        self.cutoff=cutoff
        self.index2=index2

        self.model_name='PLSuperExpCutoff'
        self.norm_par='Prefactor'

    def build_spectrum(self,prefactor_free=True,index1_free=True,cutoff_free=True,index2_free=False):
        self.parameters=[{'free':int(prefactor_free),'name':'Prefactor','min':0,
                     'max':1e6,'scale':self.prefactor_scale,'value':self.prefactor_value}]

        self.parameters.append({'free':int(index1_free),'name':'Index1','min':0,
                     'max':10,'scale':-1,'value':self.index1})

        self.parameters.append({'free':0,'name':'Scale','min':min(100,self.scale_energy*0.9),
                     'max':max(500000,self.scale_energy*1.5),'scale':1,'value':self.scale_energy})

        self.parameters.append({'free':int(cutoff_free),'name':'Cutoff','min':min(100,self.cutoff*0.9),
                     'max':max(500000,self.cutoff*1.5),'scale':1,'value':self.cutoff})

        self.parameters.append({'free':int(index2_free),'name':'Index2','min':0,'max':max(5,self.index2*1.5),
                     'value':self.index2})

        self.build(self.model_name,self.parameters)

#PLSuperExpCutoff2 class
class PLSuperExpCutoff2(Spectrum):
    def __init__(self,prefactor=1e-12,index1=2,scale_energy=1000,expfactor=1,index2=1):
        Spectrum.__init__()

        #some checks on the inputs
        if prefactor<0 or scale_energy<=0 or index2<0:
            raise ValueError("Input parameters invalid, at least one of prefactor, scale_energy, or index2 is negative or zero.")
        
        self.prefactor_scale=10**np.floor(np.log10(prefactor))
        self.prefactor=prefactor/self.prefactor_scale
        self.index1=abs(index1)
        self.scale_energy=scale_energy
        self.expfactor=expfactor/1e-3
        self.index2=index2

        self.model_name='PLSuperExpCutoff2'
        self.norm_par='Prefactor'

    def build_spectrum(self,prefactor_free=True,index1_free=True,expfactor_free=True,index2_free=False):
        self.parameters=[{'free':int(prefactor_free),'name':'Prefactor','min':0,
                     'max':1e6,'scale':self.prefactor_scale,'value':self.prefactor_value}]

        self.parameters.append({'free':int(index1_free),'name':'Index1','min':0,
                     'max':10,'scale':-1,'value':self.index1})

        self.parameters.append({'free':0,'name':'Scale','min':min(100,self.scale_energy*0.9),
                     'max':max(500000,self.scale_energy*1.5),'scale':1,'value':self.scale_energy})

        self.parameters.append({'free':int(expfactor_free),'name':'Expfactor','min':-1,
                     'max':max(100,self.expfactor*1.5),'scale':1e-3,'value':self.expfactor})

        self.parameters.append({'free':int(index2_free),'name':'Index2','min':0,'max':max(5,self.index2*1.5),
                     'value':self.index2})

        self.build(self.model_name,self.parameters)

#PLSuperExpCutoff4 class
class PLSuperExpCutoff4(Spectrum):
    def __init__(self,prefactor=1e-12,indexs=2,scale_energy=1000,expfactors=1,index2=1):
        Spectrum.__init__()

        #some checks on the inputs
        if prefactor<0 or scale_energy<=0 or index2<0:
            raise ValueError("Input parameters invalid, at least one of prefactor, scale_energy, or index2 is negative or zero.")
        
        self.prefactor_scale=10**np.floor(np.log10(prefactor))
        self.prefactor=prefactor/self.prefactor_scale
        self.indexs=abs(indexs)
        self.scale_energy=scale_energy
        self.expfactors=expfactors/1e-1
        self.index2=index2

        self.model_name='PLSuperExpCutoff4'
        self.norm_par='Prefactor'

    def build_spectrum(self,prefactor_free=True,indexs_free=True,expfactors_free=True,index2_free=False):
        self.parameters=[{'free':int(prefactor_free),'name':'Prefactor','min':0,
                     'max':1e6,'scale':self.prefactor_scale,'value':self.prefactor_value}]

        self.parameters.append({'free':int(indexs_free),'name':'IndexS','min':0,
                     'max':10,'scale':-1,'value':self.indexs})

        self.parameters.append({'free':0,'name':'Scale','min':min(100,self.scale_energy*0.9),
                     'max':max(500000,self.scale_energy*1.5),'scale':1,'value':self.scale_energy})

        self.parameters.append({'free':int(expfactors_free),'name':'ExpfactorS','min':-10,
                     'max':max(100,self.expfactors*1.5),'scale':1e-1,'value':self.expfactors})

        self.parameters.append({'free':int(index2_free),'name':'Index2','min':0,'max':max(5,self.index2*1.5),
                     'value':self.index2})

        self.build(self.model_name,self.parameters)

#LogParabola class
class LogParabola(Spectrum):
    def __init__(self,norm=1e-9,alpha=1,beta=2,eb=300):
        Spectrum.__init__()

        #some checks on the inputs
        if norm<0 or alpha<0 or eb<=0:
            raise ValueError("Input parameters invalid, at least one of norm, alpha, or eb is negative or zero.")
        
        self.norm_scale=10**np.floor(np.log10(norm))
        self.norm=norm/self.norm_scale
        self.alpha=alpha
        self.beta=beta
        self.eb=eb

        self.model_name='LogParabola'
        self.norm_par='norm'

    def build_spectrum(self,norm_free=True,alpha_free=True,beta_free=True,eb_free=False):
        self.parameters=[{'free':int(norm_free),'name':'norm','min':0,
                     'max':1e6,'scale':self.norm_scale,'value':self.norm}]

        self.parameters.append({'free':int(alpha_free),'name':'alpha','min':0,
                     'max':max(5,self.alpha*1.5),'scale':1,'value':self.alpha})

        self.parameters.append({'free':int(beta_free),'name':'beata','min':min(-5,self.beta*1.5),
                     'max':max(10,self.beta*1.5),'scale':1,'value':self.beta})

        self.parameters.append({'free':int(eb_free),'name':'Eb','min':min(10,self.eb*0.9),
                           'max':max(500000,self.eb*1.5),'scale':1,'value':self.eb})

        self.build(self.model_name,self.parameters)

#class for FileFunction
class FileFunction(Spectrum):
    def __init__(self,spectrum_file,normalization=1,apply_edisp='false'):
        Spectrum.__init__()
        
        #come checks on the inputs
        if normalization<0:
            raise ValueError(f"Input value of normalization = {normalization} is invalid, must be >0.")

        self.normalization=normalization
        self.spectrum_file=spectrum_file
        self.apply_edisp=apply_edisp

        self.model_name='FileFunction'
        self.norm_par='Normalization'

    def build_spectrum(self,normalization_free=True):
        self.parameters=[{'free':int(normalization_free),'name':'Normalization',
                    'min':min(0.01,self.nomralization*0.9),'max':max(10,self.normalization*1.5),
                    'scale':1,'value':self.normalization}]

        self.build(self.model_name,self.parameters)

        self.spectrum.setAttribute('apply_edisp',self.apply_edisp)
        self.spectrum.setAttribute('file',self.spectrum_file)

###########################################################################################
#now we'll deal with the spatial models
#not sure how to make a base class for these
###########################################################################################

#class for point sources
class SkyDir:
    def __init__(self,RA,DEC):
        #do some sanity checks
        if abs(RA)>360:
            raise ValueError(f'Input RA value of {RA} is invalid, must be between -360 and +360.')
        if abs(DEC)>90:
            raise ValueError(f'Input DEC value of {DEC} is invalid, must be between -90 and 90.')
        
        self.RA=RA
        self.DEC=DEC
        self.build()

    def build(self):
        self.document=minidom.getDOMImplementation().createDocument(None,None,None)

        self.spatial=self.document.createElement('spatialModel')
        self.spatial.setAttribute('type','SkyDirFunction')

        self.spatial.appendChild(parameter_element(free='0',name='RA',maximum='360.0',
                minimum='-360.0',scale='1.0',value=f'{self.RA:.4f}'))
        self.spatial.appendChild(parameter_element(free='0',name='DEC',maximum='90.0',
                minimum='-90.0',scale='1.0',value=f'{self.DEC:.4f}'))

#class for the RadialDisk and RadialGaussian models
class Radial:
    def __init__(self,RA,DEC,extent,extent_name):
        #do some sanity checks
        if abs(RA)>360:
            raise ValueError(f'Input RA value of {RA} is invalid, must be between -360 and +360.')
        if abs(DEC)>90:
            raise ValueError(f'Input DEC value of {DEC} is invalid, must be between -90 and 90.')
        if extent<=0:
            raise ValueError(f'Input extent value of {extent} is invalid, must be >0.')
        if extent_name not in ['Radius','Sigma']:
            raise ValueError(f'Input name for extent parameter is invalid, must be "Radius" or "Sigma".')
        
        self.RA=RA
        self.DEC=DEC
        self.extent=extent
        self.extent_name=extent_name
        self.model_name='RadialDisk' if extent_name=='Radius' else 'RadialGaussian'

        self.build()

    def build(self):
        self.document=minidom.getDOMImplementation().createDocument(None,None,None)

        self.spatial=self.document.createElement('spatialModel')
        self.spatial.setAttribute('type',self.model_name)

        self.spatial.appendChild(parameter_element(free='0',name=self.extent_name,
                maximum=max(10,self.extent*1.5),minimum='0',scale='1',value=self.extent))
        self.spatial.appendChild(parameter_element(free='0',name='RA',maximum='360.0',
                minimum='-360.0',scale='1.0',value=f'{self.RA:.4f}'))
        self.spatial.appendChild(parameter_element(free='0',name='DEC',maximum='90.0',
                minimum='-90.0',scale='1.0',value=f'{self.DEC:.4f}'))

#class for extended sources with spatial template file
class SpatialMap:
    def __init__(self,spatial_file):
        self.spatial_file=spatial_file

        self.build()

    def build(self):
        self.document=minidom.getDOMImplementation().createDocument(None,None,None)

        self.spatial=self.document.createElement('spatialModel')
        self.spatial.setAttribute('type','SpatialMap')
        self.spatial.setAttribute('file',self.spatial_file)
        self.spatial.setAttribute('map_based_integral','true')

        self.spatial.appendChild(parameter_element(free="0",name="Prefactor",maximum="1000",
                minimum="0.001",scale="1",value="1"))

#class for constant value spatial type
class ConstantValue:
    def __init__(self,value=1):
        #check on inputs
        if value<0:
            raise ValueError(f'Input value of {value} is invalid, must be >=0.')

        self.value=value

        self.build()

    def build(self):
        self.document=minidom.getDOMImplementation().createDocument(None,None,None)

        self.spatial=self.document.createElement('spatialModel')
        self.spatial.setAttribute('type','ConstantValue')

        self.spatial.appendChild(parameter_element(free="0",name="Value",maximum=max(10,self.value*1.5),
                    minimum='0',scale='1',value=self.value))

#class for mapcubefunction spatial model
class MapCubeFunction
    def __init__(self,spatial_file,normalization=1):
        #check on input value
        if normalization<=0:
            raise ValueError(f'Input value for normalizaton of {normalization} is invalid, must be >0.')
        
        self.spatial_file=spatial_file
        self.normalization=normalization

        self.build()

    def build(self):
        self.document=minidom.getDOMImplementation().createDocument(None,None,None)

        self.spatial=self.document.createElement('spatialModel')
        self.spatial.setAttribute('type','MapCubeFunction')
        self.spatial.setAttribute('file',self.spatial_file)

        self.spatial.appendChild(parameter_element(free="0",name="Normalization",
                maximum=max(1000,self.normalization*1.5),
                minimum=min(0.001,self.normalization*0.9),scale=1,value=self.normalization))



        
