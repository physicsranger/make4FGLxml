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
#class for the spectral models, will need a function for each model
#currently, focus mainly on those models we'll run into in the catalog
#only need to pass the spectrum attribute to the appendChild method of a source XML element
###########################################################################################
class Spectrum:
    def __init__(self,model,**kwargs):
        self.document=minidom.getDOMImplementation().createDocument(None,None,None)
        self.spectrum=self.document.createElement('spectrum')

        self.get_function_dictionary()

        self.parameters=self.functions[model]

        self.model_name=model

    def build(self):
        self.spectrum.setAttribute('type',self.model_name)
        for parameter in self.parameters:
            self.spectrum.appendChild(parameter_element(parameter.get('free'),
                parameter.get('name'),parameter.get('max'),parameter.get('min'),
                parameter.get('scale'),parameter.get('value')))

    def get_function_dictionary(self):
        self.functions={'PowerLaw':self.PowerLaw,
                   'PowerLaw2':self.PowerLaw2,
                   'PLSuperExpCutoff':self.PLSuperExpCutoff,
                   'PLSuperExpCutoff2':self.PLSuperExpCutoff2,
                   'PLSuperExpCutoff4':self.PLSuperExpCutoff4,
                   'LogParabola':self.LogParabola,
                   'FileFunction':self.FileFunction}
                   

    def PowerLaw(self,Prefactor=1e-12,Scale=1000,Index2,Prefactor_free=True,Index_free=True,model=None):
        if Prefactor<0 or Scale<=0:
            raise ValueError("Input parameters invalid, at least one of Prefactor or Scale is negative or zero.")
        
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=prefactor/self.Prefactor_scale
        self.Index=abs(Index)
        self.Scale=scale

        self.model_name='PowerLaw'
        self.norm_par='Prefactor'

        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','min':0,
                     'max':1e6,'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(Index_free),'name':'Index','min':0,
                     'max':10,'scale':-1,'value':self.Index})

        self.parameters.append({'free':0,'name':'Scale','min':min(100,self.Scale*0.9),
                           'max':max(500000,self.Scale*1.5),'scale':1,'value':self.Scale})

        self.build(self.model_name,self.parameters)

    def PowerLaw2(self,Integral=1e-8,Index=2,LowerLimit=100,UpperLimit=300,
                  Integral_free=True,Index_free=True,model=None):
        #some checks on the inputs
        if Integral<0 or LowerLimit<=0 or UpperLimit<=0:
            raise ValueError("Input parameters invalid, at least one of Integral, LowerLimit, or UpperLimit is negative or zero.")
        
        if LowerLimit>=UpperLimit:
            raise ValueError(f"Requested LowerLimit of {LowerLimit:,} is >= UpperLimit of {UpperLimit:,}.")
        
        self.Integral_scale=10**np.floor(np.log10(Integral))
        self.Integral=Integral/self.Integral_scale
        self.index=abs(index)
        self.LowerLimit=LowerLimit
        self.UpperLimit=UpperLimit

        self.model_name='PowerLaw2'
        self.norm_par='Integral'

        self.parameters=[{'free':int(Integral_free),'name':'Integral','min':0,
                     'max':1e6,'scale':self.Integral_scale,'value':self.Integral}]

        self.parameters.append({'free':int(Index_free),'name':'Index','min':0,
                     'max':10,'scale':-1,'value':self.Index})

        self.parameters.append({'free':0,'name':'LowerLimit','min':10,
                           'max':max(500000,self.LowerLimit*1.5),
                           'scale':1,'value':self.LowerLimit})

        self.parameters.append({'free':0,'name':'UpperLimit','min':10,
                           'max':max(500000,self.UpperLimit*1.5),
                           'scale':1,'value':self.UpperLimit})

        self.build(self.model_name,self.parameters)

    def PLSuperExpCutoff(self,Prefactor=1e-12,Index1=2,Scale=1000,Cutoff=1000,Index2=1,
                         Prefactor_free=True,Index1_free=True,Cutoff_free=True,Index2_free=False,model=None):
        #some checks on the inputs
        if Prefactor<0 or Scale<=0 or Cutoff<=0 or Index2<0:
            raise ValueError("Input parameters invalid, at least one of Prefactor, Scale, Cutoff, or Index2 is negative or zero.")
        
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=Prefactor/self.Prefactor_scale
        self.Index1=abs(Index1)
        self.Scale=Scale
        self.Cutoff=Cutoff
        self.Index2=Index2

        self.model_name='PLSuperExpCutoff'
        self.norm_par='Prefactor'

        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','min':0,
                     'max':1e6,'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(index1_free),'name':'Index1','min':0,
                     'max':10,'scale':-1,'value':self.Index1})

        self.parameters.append({'free':0,'name':'Scale','min':min(100,self.Scale*0.9),
                     'max':max(500000,self.Scale*1.5),'scale':1,'value':self.Scale})

        self.parameters.append({'free':int(Cutoff_free),'name':'Cutoff','min':min(100,self.Cutoff*0.9),
                     'max':max(500000,self.Cutoff*1.5),'scale':1,'value':self.Cutoff})

        self.parameters.append({'free':int(Index2_free),'name':'Index2','min':0,'max':max(5,self.Index2*1.5),
                     'value':self.Index2})

        self.build(self.model_name,self.parameters)

    def PLSuperExpCutoff2(self,Prefactor=1e-12,Index1=2,Scale=1000,Expfactor=1,Index2=1,
                          Prefactor_free=True,Index1_free=True,Expfactor_free=True,Index2_free=False,model=None):
        #some checks on the inputs
        if Prefactor<0 or Scale<=0 or Index2<0:
            raise ValueError("Input parameters invalid, at least one of Prefactor, Scale, or Index2 is negative or zero.")
        
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=Prefactor/self.Prefactor_scale
        self.Index1=abs(Index1)
        self.Scale=Scale
        self.Expfactor=Expfactor/1e-3
        self.Index2=Index2

        self.model_name='PLSuperExpCutoff2'
        self.norm_par='Prefactor'

        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','min':0,
                     'max':1e6,'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(Index1_free),'name':'Index1','min':0,
                     'max':10,'scale':-1,'value':self.Index1})

        self.parameters.append({'free':0,'name':'Scale','min':min(100,self.Scale*0.9),
                     'max':max(500000,self.Scale*1.5),'scale':1,'value':self.Scale})

        self.parameters.append({'free':int(Expfactor_free),'name':'Expfactor','min':-1,
                     'max':max(100,self.Expfactor*1.5),'scale':1e-3,'value':self.Expfactor})

        self.parameters.append({'free':int(Index2_free),'name':'Index2','min':0,'max':max(5,self.Index2*1.5),
                     'value':self.Index2})

        self.build(self.model_name,self.parameters)

    def PLSuperExpCutoff4(self,refactor=1e-12,IndexS=2,Scale=1000,ExpfactorS=1,Index2=1,
                          Prefactor_free=True,IndexS_free=True,ExpfactorS_free=True,Index2_free=False,model=None):
        #some checks on the inputs
        if Prefactor<0 or Scale<=0 or Index2<0:
            raise ValueError("Input parameters invalid, at least one of Prefactor, Scale, or Index2 is negative or zero.")
        
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=Prefactor/self.Prefactor_scale
        self.IndexS=abs(IndexS)
        self.Scale=Scale
        self.ExpfactorS=ExpfactorS/1e-1
        self.Index2=Index2

        self.model_name='PLSuperExpCutoff4'
        self.norm_par='Prefactor'

        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','min':0,
                     'max':1e6,'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(IndexS_free),'name':'IndexS','min':0,
                     'max':10,'scale':-1,'value':self.IndexS})

        self.parameters.append({'free':0,'name':'Scale','min':min(100,self.Scale*0.9),
                     'max':max(500000,self.Scale*1.5),'scale':1,'value':self.Scale})

        self.parameters.append({'free':int(ExpfactorS_free),'name':'ExpfactorS','min':-10,
                     'max':max(100,self.ExpfactorS*1.5),'scale':1e-1,'value':self.ExpfactorS})

        self.parameters.append({'free':int(Index2_free),'name':'Index2','min':0,'max':max(5,self.Index2*1.5),
                     'value':self.Index2})

        self.build(self.model_name,self.parameters)

    def LogParabola(self,norm=1e-9,alpha=1,beta=2,Eb=300,
                    norm_free=True,alpha_free=True,beta_free=True,Eb_free=False,model=None):
        #some checks on the inputs
        if norm<0 or alpha<0 or Eb<=0:
            raise ValueError("Input parameters invalid, at least one of norm, alpha, or Eb is negative or zero.")
        
        self.norm_scale=10**np.floor(np.log10(norm))
        self.norm=norm/self.norm_scale
        self.alpha=alpha
        self.beta=beta
        self.Eb=Eb

        self.model_name='LogParabola'
        self.norm_par='norm'

        self.parameters=[{'free':int(norm_free),'name':'norm','min':0,
                     'max':1e6,'scale':self.norm_scale,'value':self.norm}]

        self.parameters.append({'free':int(alpha_free),'name':'alpha','min':0,
                     'max':max(5,self.alpha*1.5),'scale':1,'value':self.alpha})

        self.parameters.append({'free':int(beta_free),'name':'beata','min':min(-5,self.beta*1.5),
                     'max':max(10,self.beta*1.5),'scale':1,'value':self.beta})

        self.parameters.append({'free':int(Eb_free),'name':'Eb','min':min(10,self.Eb*0.9),
                           'max':max(500000,self.Eb*1.5),'scale':1,'value':self.Eb})

        self.build(self.model_name,self.parameters)

    def FileFunction(self,spectrum_file,normalization=1,apply_edisp='false',
                     normalization_free=True,model=None):
        #come checks on the inputs
        if normalization<0:
            raise ValueError(f"Input value of normalization = {normalization} is invalid, must be >0.")

        self.normalization=normalization
        self.spectrum_file=spectrum_file
        self.apply_edisp=apply_edisp

        self.model_name='FileFunction'
        self.norm_par='Normalization'

        self.parameters=[{'free':int(normalization_free),'name':'Normalization',
                    'min':min(0.01,self.nomralization*0.9),'max':max(10,self.normalization*1.5),
                    'scale':1,'value':self.normalization}]

        self.build(self.model_name,self.parameters)

        self.spectrum.setAttribute('apply_edisp',self.apply_edisp)
        self.spectrum.setAttribute('file',self.spectrum_file)

###########################################################################################
#class for the spatial models
#will only need to add the spatial class attribute/object to a larger XML document
###########################################################################################
class Spatial:
    def __init__(self,spatial_model,**kwargs):
        self.document=minidom.getDOMImplementation().createDocument(None,None,None)
        self.spatial=self.document.createElement('spatialModel')

        self.get_function_dictionary()

        self.functions[spatial_model](**kwargs)

    def get_function_dictionary(self):
        self.functions={'SkyDir':self.SkyDir,
                   'RadialDisk':self.RadialDisk,
                   'RadialGaussian':self.RadialGaussian,
                   'SpatialMap':self.SpatialMap,
                   'ConstantValue':self.ConstantValue,
                   'MapCubeFunction':self.MapCubeFunction}

    def SkyDir(self,RA,DEC,spatial_model=None):
        #do some sanity checks
        if abs(RA)>360:
            raise ValueError(f'Input RA value of {RA} is invalid, must be between -360 and +360.')
        if abs(DEC)>90:
            raise ValueError(f'Input DEC value of {DEC} is invalid, must be between -90 and 90.')
        
        self.RA=RA
        self.DEC=DEC
        self.build()

        self.spatial.setAttribute('type','SkyDirFunction')

        self.spatial.appendChild(parameter_element(free='0',name='RA',maximum='360.0',
                minimum='-360.0',scale='1.0',value=f'{self.RA:.4f}'))
        self.spatial.appendChild(parameter_element(free='0',name='DEC',maximum='90.0',
                minimum='-90.0',scale='1.0',value=f'{self.DEC:.4f}'))

    
    def RadialDisk(self,RA,DEC,Radius,spatial_model=None):
        #do some sanity checks
        if abs(RA)>360:
            raise ValueError(f'Input RA value of {RA} is invalid, must be between -360 and +360.')
        if abs(DEC)>90:
            raise ValueError(f'Input DEC value of {DEC} is invalid, must be between -90 and 90.')
        if Radius<=0:
            raise ValueError(f'Input Radius value of {Radius} is invalid, must be >0.')
        
        self.RA=RA
        self.DEC=DEC
        self.Radius=Radius

        self.spatial.setAttribute('type','RadialDisk')

        self.spatial.appendChild(parameter_element(free='0',name='Radius',
                maximum=max(10,self.extent*1.5),minimum='0',scale='1',value=self.Radius))
        self.spatial.appendChild(parameter_element(free='0',name='RA',maximum='360.0',
                minimum='-360.0',scale='1.0',value=f'{self.RA:.4f}'))
        self.spatial.appendChild(parameter_element(free='0',name='DEC',maximum='90.0',
                minimum='-90.0',scale='1.0',value=f'{self.DEC:.4f}'))

    def RadialGaussian(self,RA,DEC,Sigma,spatial_model=None):
        #do some sanity checks
        if abs(RA)>360:
            raise ValueError(f'Input RA value of {RA} is invalid, must be between -360 and +360.')
        if abs(DEC)>90:
            raise ValueError(f'Input DEC value of {DEC} is invalid, must be between -90 and 90.')
        if Radius<=0:
            raise ValueError(f'Input Sigma value of {Sigma} is invalid, must be >0.')
        
        self.RA=RA
        self.DEC=DEC
        self.Sigma=Sigma

        self.spatial.setAttribute('type','RadialGaussian')

        self.spatial.appendChild(parameter_element(free='0',name='Sigma',
                maximum=max(10,self.extent*1.5),minimum='0',scale='1',value=self.Sigma))
        self.spatial.appendChild(parameter_element(free='0',name='RA',maximum='360.0',
                minimum='-360.0',scale='1.0',value=f'{self.RA:.4f}'))
        self.spatial.appendChild(parameter_element(free='0',name='DEC',maximum='90.0',
                minimum='-90.0',scale='1.0',value=f'{self.DEC:.4f}'))

    def SpatialMap(self,spatial_file,spatial_model=None):
        self.spatial_file=spatial_file

        self.spatial.setAttribute('type','SpatialMap')
        self.spatial.setAttribute('file',self.spatial_file)
        self.spatial.setAttribute('map_based_integral','true')

        self.spatial.appendChild(parameter_element(free="0",name="Prefactor",maximum="1000",
                minimum="0.001",scale="1",value="1"))

    def ConstantValue(self,value=1,spatial_model=None):
        #check on inputs
        if value<0:
            raise ValueError(f'Input value of {value} is invalid, must be >=0.')

        self.value=value

        self.spatial.setAttribute('type','ConstantValue')

        self.spatial.appendChild(parameter_element(free="0",name="Value",maximum=max(10,self.value*1.5),
                    minimum='0',scale='1',value=self.value))

    def MapCubeFunction(self,spatial_file,normalization=1,spatial_model=None):
        #check on input value
        if normalization<=0:
            raise ValueError(f'Input value for normalizaton of {normalization} is invalid, must be >0.')
        
        self.spatial_file=spatial_file
        self.normalization=normalization

        self.spatial.setAttribute('type','MapCubeFunction')
        self.spatial.setAttribute('file',self.spatial_file)

        self.spatial.appendChild(parameter_element(free="0",name="Normalization",
                maximum=max(1000,self.normalization*1.5),
                minimum=min(0.001,self.normalization*0.9),scale=1,value=self.normalization))



        
