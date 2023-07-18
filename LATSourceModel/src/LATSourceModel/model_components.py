from xml.dom import minidom

from LATSourceModel.utilities import parameter_element

import os

import numpy as np



###########################################################################################
#class for the spectral models
###########################################################################################
class Spectrum:
    '''
    A class to create a spectrum XML document element to be added to a source
    document element in a Fermi LAT XML region model

    ...

    Attributes
    ----------
    document : XML document
        document instance containing the spectrum document element
    functions : dict
        dictionary with keys as spectral function names and
        values pointing to the corresponding functions
    model_name : str
        name of spectrum function, e.g., PowerLaw
    norm_par : str
        name of spectral model normalization parameter
    parameters : list
        a list of dictionaries, with each dictionary containing
        the necessary keys and values to create a parameter_element
        (see build_model.utilities for more information)
    spectrum : XML document element
        document element corresponding to the spectrum model

    Methods
    -------
    build()
        build the spectrum document element, adding all of the
        parameters in the parameters attribute
    get_function_dictionary()
        create the functions attribute as a dictionary of spectral model
        functions
    
    <for descriptions of all spectral models referenced in the methods below
    see https://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/source_models.html>
    PowerLaw(Prefactor=1e-12,Scale=1000,Index=2,Prefactor_free=True,
      Index_free=True,Scale_free=False,model=None)
        create parameter dictionaries corresponding to the PowerLaw spectral model
    LogParabola(self,norm=1e-9,alpha=1,beta=2,Eb=300,norm_free=True,
      alpha_free=True,beta_free=True,Eb_free=False,model=None)
        create parameter dictionaries corresponding to LogParabola spectral model
    PLSuperExpCutoff4(self,Prefactor=1e-12,IndexS=2,Scale=1000,ExpfactorS=1,
      Index2=1,Prefactor_free=True,IndexS_free=True,ExpfactorS_free=True,
      Index2_free=False,Scale_free=False,model=None)
        create parameter dictionaries corresponding to PLSuperExpCutoff4 spectral model
    PLSuperExpCutoff2(self,Prefactor=1e-12,Index1=2,Scale=1000,Expfactor=1,
      Index2=1,Prefactor_free=True,Index1_free=True,Expfactor_free=True,
      Index2_free=False,Scale_free=False,model=None)
        create parameter dictionaries corresponding to PLSuperExpCutoff2 spectral model
    PowerLaw2(self,Integral=1e-8,Index=2,LowerLimit=100,UpperLimit=300000,
      Integral_free=True,Index_free=True,model=None)
        create parameter dictionaries corresponding to PowerLaw2 spectral model

    BandFunction(self,norm=1e-9,alpha=1.8,beta=2.5,Ep=0.1,norm_free=True,
      alpha_free=True,beta_free=True,Ep_free=True,model=None)
        create parameter dictionaries corresponding to BandFunction spectral model
    BPLExpCutoff(self,Prefactor=1e-12,Index1=1.8,Index2=2.3,BreakValue=1000,
      Eabs=10,P1=100,Prefactor_free=True,Index1_free=True,Index2_free=False,
      BreakValue_free=True,Eabs_free=True,P1_free=False,model=None)
        create parameter dictionaries corresponding to BPLExpCutoff spectral model
    BrokenPowerLaw(self,Prefactor=1e-12,Index1=1.8,BreakValue=1000,
      Index2=2.3,Prefactor_free=True,Index1_free=True,BreakValue_free=False,
      Index2_free=True,model=None)
        create parameter dictionaries corresponding to the BrokenPowerLaw spectral model
    BrokenPowerLaw2(self,Integral=1e-8,Index1=1.8,BreakValue=1000,Index2=2.3,
      LowerLimit=100,UpperLimit=300000,Integral_free=True,Index1_free=True,
      BreakValue_free=False,Index2_free=True,model=None)
        create parameter dictionaries corresponding to BrokenPowerLaw2 spectral model
    ConstantValue(self,Value=1,Value_free=True,model=None)
        create parameter dictionaries corresponding to ConstantValue spectral model
    DMFitFunction(self,
      spectrum_file=os.path.join('$(FERMI_DIR)','refdata','fermi','Likelihood','gammamc_dif.dat'),
      norm=5e20,sigmav=3e-26,mass=10,bratio=1,channel0=4,channel1=1,
      norm_free=True,sigmav_free=True,mass_free=True,bratio_free=False,
      channel0_free=False,channel1_free=False,model=None)
        create parameter dictionaries corresponding to DMFitFunction spectral model
    ExpCutoff(self,Prefactor=1e-12,Index=2,Scale=1000,Ebreak=1000,P1=100,P2=0,P3=0,
      Prefactor_free=True,Index_free=True,Ebreak_free=True,P1_free=False,
      P2_free=False,P3_free=False,Scale_free=False,model=None)
        create parameter dictionaries corresponding to ExpCutoff spectral model
    FileFunction(self,spectrum_file,Normalization=1,apply_edisp='false',
      Normalization_free=True,model=None)
        create parameter dictionaries corresponding to FileFunction spectral model
    Gaussian(self,Prefactor=1e-9,Mean=7e4,Sigma=1e3,Prefactor_free=True,
      Mean_free=True,Sigma_free=True,model=None)
        create parameter dictionaries corresponding to Gaussian spectral model
    PLSuperExpCutoff(self,Prefactor=1e-12,Index1=2,Scale=1000,Cutoff=1000,
      Index2=1,Prefactor_free=True,Index1_free=True,Cutoff_free=True,
      Index2_free=False,Scale_free=False,model=None)
        create parameter dictionaries correspnding to PLSuperExpCutoff spectral model
    PLSuperExpCutoff3(self,Prefactor=1e-12,IndexS=2,Scale=1000,Expfactor2=2,
      Index2=1,Prefactor_free=True,IndexS_free=True,Expfactor2_free=True,
      Index2_free=False,Scale_free=False,model=None)
        create parameter dictionaries corresponding to PLSuperExpCutoff3 spectral model
    SmoothBrokenPowerLaw(self,Prefactor=1e-12,Index1=1.8,Scale=100,BreakValue=1000,
      Index2=2.3,Beta=0.2,Prefactor_free=True,Index1_free=True,Scale_free=False,
      BreakValue_free=True,Index2_free=True,Beta_free=True,model=None)
        create parameter dictionaries corresponding to SmoothBrokenPowerLaw spectral model
    
    Raises
    ------
    ValueError
        if requested spectral model in not in the function dictionary,
        raise an error
    '''
    
    def __init__(self,model,**kwargs):
        '''
        Parameters
        ----------
        model - str
            name of spectral model
        **kwargs - various
            named arguments matching the arguments of function corresponding
            to model, typically passed in as a dictionary
        '''
        
        self.document=minidom.getDOMImplementation().createDocument(None,None,None)
        self.spectrum=self.document.createElement('spectrum')

        self.get_function_dictionary()

         #do some checks on the input
        if model not in self.functions.keys():
            raise ValueError('{model = } not implemented')

        self.parameters=self.functions[model](**kwargs)

    def build(self):
        '''
        construct XML spectrum document element
        '''

        #set basic spectrum attributes
        self.spectrum.setAttribute('type',self.model_name)
        self.spectrum.setAttribute('normPar',self.norm_par)

        #cycle through parameter dictionaries
        #adding elements to spectrum document element
        for parameter in self.parameters:
            self.spectrum.appendChild(parameter_element(**parameter))

    def get_function_dictionary(self):
        '''
        construct functions attribute as a dictionary
        of available spectral model functions
        '''
        
        self.functions={'PowerLaw':self.PowerLaw,
                   'BrokenPowerLaw':self.BrokenPowerLaw,
                   'PowerLaw2':self.PowerLaw2,
                   'BrokenPowerLaw2':self.BrokenPowerLaw2,
                   'SmoothBrokenPowerLaw':self.SmoothBrokenPowerLaw,
                   'ExpCutoff':self.ExpCutoff,
                   'BPLExpCutoff':self.BPLExpCutoff,
                   'Gaussian':self.Gaussian,
                   'ConstantValue':self.ConstantValue,
                   'BandFunction':self.BandFunction,
                   'PLSuperExpCutoff':self.PLSuperExpCutoff,
                   'PLSuperExpCutoff2':self.PLSuperExpCutoff2,
                   'PLSuperExpCutoff3':self.PLSuperExpCutoff3,
                   'PLSuperExpCutoff4':self.PLSuperExpCutoff4,
                   'LogParabola':self.LogParabola,
                   'DMFitFunction':self.DMFitFunction,
                   'FileFunction':self.FileFunction}
                   

    def PowerLaw(self,Prefactor=1e-12,Scale=1000,Index=2,Prefactor_free=True,Index_free=True,Scale_free=False,model=None):
        '''
        create parameter dictionaries for PowerLaw spectral model and then
        build spectrum document element

        Parameters
        ----------
        Prefactor - float
            value of Prefactor parameter, normalization in MeV**-1 cm**-2 s**-1
        Scale - float
            value of Scale parameter, scale energy in MeV
        Index - float
            value of Index parameter, spectral index using E**-(Index) convention
        Prefactor_free - bool
            flag to set the Prefactor parameter free
        Index_free - bool
            flag to set the Index parameter free
        Scale_free - bool
            flag to set the Scale parameter free, NOTE: it is not
            recommended to set this to True
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if Prefactor and Scale inputs are invalid choices, raise an error
        '''

        #some checks on input values
        if Prefactor<0 or Scale<=0:
            raise ValueError("Input parameters invalid, at least one of Prefactor or Scale is negative or zero.")

        #massage the inputs a little and save
        #parameters as attributes to be accesible later
        #if ever desired
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=Prefactor/self.Prefactor_scale
        self.Index=abs(Index)
        self.Scale=Scale

        #set descriptive attributes
        self.model_name='PowerLaw'
        self.norm_par='Prefactor'

        #make a dictionary for each parameter
        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','minimum':0,
                     'maximum':1e6,'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(Index_free),'name':'Index','minimum':0,
                     'maximum':10,'scale':-1,'value':self.Index})

        self.parameters.append({'free':int(Scale_free),'name':'Scale','minimum':min(100,self.Scale*0.9),
                           'maximum':max(500000,self.Scale*1.5),'scale':1,'value':self.Scale})

        #now build the spectrum document element
        self.build()

    def BrokenPowerLaw(self,Prefactor=1e-12,Index1=1.8,BreakValue=1000,Index2=2.3,
                       Prefactor_free=True,Index1_free=True,BreakValue_free=False,
                       Index2_free=True,model=None):
        '''
        create parameter dictionaries for BrokenPowerLaw spectral model and then
        build spectrum document element

        Parameters
        ----------
        Prefactor - float
            value of Prefactor parameter, normalization in MeV**-1 cm**-2 s**-1
        Index1 - float
            value of Index1 parameter, low-energy spectral index using E**-(Index1) convention
        Breakvalue - float
            value of BreakValue parameter, energy in MeV of spectral break
        Index2 - float
            value of Index2 parameter, high-energy spectral index using E**-(Index2) convention
        Prefactor_free - bool
            flag to set the Prefactor parameter free
        Index1_free - bool
            flag to set the Index1 parameter free
        BreakValue_free - bool
            flag to set the BreakValue parameter free
        Index2_free - bool
            flag to set the Index2 parameter free
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        If Prefactor or BreakValue inputs are invalid, raise an error
        '''

        #do some sanity checks on the input values
        if Prefactor<0 or BreakValue<=0:
            raise ValueError("Input parameters invalid, at least one of Prefactor or BreakValue is negative or zero.")

        #massage the input values a bit
        #and make accessible for later
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=Prefactor/self.Prefactor_scale
        self.Index1=abs(Index1)
        self.BreakValue=BreakValue
        self.Index2=abs(Index2)

        #set descriptive attributes
        self.model_name='BrokenPowerLaw'
        self.norm_par='Prefactor'

        #make a dictionary for each parameter
        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','minimum':0,
                     'maximum':1e6,'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(Index1_free),'name':'Index1','minimum':0,
                     'maximum':10,'scale':-1,'value':self.Index1})

        self.parameters.append({'free':int(BreakValue_free),'name':'BreakValue',
                     'minimum':min(100,self.BreakValue*0.9),
                     'maximum':max(500000,self.BreakValue*1.5),'scale':1,'value':self.BreakValue})

        self.parameters.append({'free':int(Index2_free),'name':'Index2','minimum':0,
                     'maximum':10,'scale':-1,'value':self.Index2})

        #now build the spectrum document element
        self.build()

    def PowerLaw2(self,Integral=1e-8,Index=2,LowerLimit=100,UpperLimit=300000,
                  Integral_free=True,Index_free=True,LowerLimit_free=False,
                  UpperLimit_free=False,model=None):
        '''
        create parameter dictionaries for PowerLaw2 spectral model and then
        build spectrum document element

        Parameters
        ----------
        Integral - float
            value of Integral parameter, integral flux (from LowerLimit to
            UpperLimit) in units of cm**-2 s**-1
        Index - float
            value of Index parameter, spectral index using E**-(Index) convention
        LowerLimit - float
            value of LowerLimit parameter, lower energy value, in MeV, of the
            flux integral
        UpperLimit - float
            value of UpperLimit parameter, upper energy value, in MeV, of the
            flux integral
        Integral_free - bool
            flag to set the Integral parameter free
        Index_free - bool
            flag to set the Index parameter free
        LowerLimit_free - bool
            flag to set the LowerLimit parameter free, NOTE: it is not
            recommended to set this flag to True
        UpperLimit_free - bool
            flag to set the UpperLimit parameter free, NOTE: it is not
            recommended to set this flag to True
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if Integral, LowerLimit, or UpperLimit input values are invalid,
            raise an error
        '''
        
        #some checks on the inputs
        if Integral<0 or LowerLimit<=0 or UpperLimit<=0:
            raise ValueError("Input parameters invalid, at least one of Integral, LowerLimit, or UpperLimit is negative or zero.")
        
        if LowerLimit>=UpperLimit:
            raise ValueError(f"Requested of {LowerLimit = :,} is >= {UpperLimit = :,}.")

        #massage the inputs a bit and make accessible
        #for later use, if desired
        self.Integral_scale=10**np.floor(np.log10(Integral))
        self.Integral=Integral/self.Integral_scale
        self.Index=abs(Index)
        self.LowerLimit=LowerLimit
        self.UpperLimit=UpperLimit

        #assign descriptive attributes
        self.model_name='PowerLaw2'
        self.norm_par='Integral'

        #create dictionary for each parameter
        self.parameters=[{'free':int(Integral_free),'name':'Integral','minimum':0,
                     'maximum':1e6,'scale':self.Integral_scale,'value':self.Integral}]

        self.parameters.append({'free':int(Index_free),'name':'Index','minimum':0,
                     'maximum':10,'scale':-1,'value':self.Index})

        self.parameters.append({'free':int(LowerLimit_free),'name':'LowerLimit','minimum':10,
                           'maximum':max(500000,self.LowerLimit*1.5),
                           'scale':1,'value':self.LowerLimit})

        self.parameters.append({'free':int(UpperLimit_free),'name':'UpperLimit','minimum':10,
                           'maximum':max(500000,self.UpperLimit*1.5),
                           'scale':1,'value':self.UpperLimit})

        #now buld the spectrum document element
        self.build()

    def BrokenPowerLaw2(self,Integral=1e-8,Index1=1.8,BreakValue=1000,Index2=2.3,
                       LowerLimit=100,UpperLimit=300000,Integral_free=True,Index1_free=True,
                       BreakValue_free=False,Index2_free=True,LowerLimit_free=False,
                       UpperLimit_free=False,model=None):
        '''
        create parameter dictionaries for BrokenPowerLaw2 spectral model and then
        build spectrum document element

        Parameters
        ----------
        Integral - float
            value of Integral parameter, integrated flux, from LowerLimit
            to UpperLimit, in units of cm**-2 s**-1
        Index1 - float
            value of Index1 parameter, low-energy spectral index using E**-(Index1) convention
        BreakValue - float
            value of the BreakValue parameter, energy in MeV of spectral break
        Index2 - float
            value of Index2 parameter, high-energy spectral index using E**-(Index2) convention
        LowerLimit - float
            value of LowerLimit parameter, low energy value, in MeV, of the
            flux integral
        UpperLimit - float
            value of UpperLimit parameter, upper energy value, in MeV, of the
            flux integral
        Integral_free - bool
            flag to set the Integral parameter free
        Index1_free - bool
            flag to set the Index1 parameter free
        BreakValue_free - bool
            flag to set the BreakValue parameter free
        Index2_free - bool
            flag to set the Index2 parameter free
        LowerLimit_free - bool
            flag to set the LowerLimit parameter free, NOTE: it is not
            recommended to set this flag to True
        UpperLimit_free - bool
            flag to set the UpperLimit parameter free, NOTE: it is not
            recommended to set this flag to True
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if Integral, LowerLimit, UpperLimit, or BreakValue inputs are
            invalid, raise an error
        '''
        
        #some checks on the inputs
        if Integral<0 or LowerLimit<=0 or UpperLimit<=0 or BreakValue<=0:
            raise ValueError("Input parameters invalid, at least one of Integral, LowerLimit, UpperLimit, or BreakValue is negative or zero.")
        
        if LowerLimit>=UpperLimit:
            raise ValueError(f"Requested of {LowerLimit = :,} is >= {UpperLimit = :,}.")

        #massage the inputs a little and make
        #them accessible later, if desired
        self.Integral_scale=10**np.floor(np.log10(Integral))
        self.Integral=Integral/self.Integral_scale
        self.Index1=abs(Index1)
        self.BreakValue=BreakValue
        self.Index2=abs(Index2)
        self.LowerLimit=LowerLimit
        self.UpperLimit=UpperLimit

        #set some descriptive attributes
        self.model_name='BrokenPowerLaw2'
        self.norm_par='Integral'

        #create dictionary for each parameter
        self.parameters=[{'free':int(Integral_free),'name':'Integral','minimum':0,
                     'maximum':1e6,'scale':self.Integral_scale,'value':self.Integral}]

        self.parameters.append({'free':int(Index1_free),'name':'Index1','minimum':0,
                     'maximum':10,'scale':-1,'value':self.Index1})

        self.parameters.append({'free':int(BreakValue_free),'name':'BreakValue',
                     'minimum':min(100,self.BreakValue*0.9),
                     'maximum':max(500000,self.BreakValue*1.5),'scale':1,'value':self.BreakValue})

        self.parameters.append({'free':int(Index2_free),'name':'Index2','minimum':0,
                     'maximum':10,'scale':-1,'value':self.Index2})

        self.parameters.append({'free':int(LowerLimit_free),'name':'LowerLimit','minimum':10,
                           'maximum':max(500000,self.LowerLimit*1.5),
                           'scale':1,'value':self.LowerLimit})

        self.parameters.append({'free':int(UpperLimit_free),'name':'UpperLimit','minimum':10,
                           'maximum':max(500000,self.UpperLimit*1.5),
                           'scale':1,'value':self.UpperLimit})

        #now build the spectrum document element
        self.build()

    def SmoothBrokenPowerLaw(self,Prefactor=1e-12,Index1=1.8,Scale=100,BreakValue=1000,
                             Index2=2.3,Beta=0.2,Prefactor_free=True,Index1_free=True,
                             Scale_free=False,BreakValue_free=True,Index2_free=True,
                             Beta_free=True,model=None):
        '''
        create parameter dictionaries for SmoothBrokenPowerLaw model and then
        build spectrum document element

        Parameters
        ----------
        Prefactor - float
            value of Prefactor parameter, normalization in MeV**-1 cm**-2 s**-1
        Index1 - float
            value of Index1 parameter, low-energy spectral index using E**-(Index1) convention
        Scale - float
            value of Scale parameter, scale energy in MeV
        BreakValue - float
            value of BreakValue parameter, spectral break energy in MeV
        Index2 - float
            value of Index2 parameter, high_energy spectral index using E**-(Index2) convention
        Beta - float
            value of Beta parameter
        Prefactor_free - bool
            flag to set the Prefactor parameter free
        Index1_free - bool
            flag to set the Index1 parameter free
        Scale_free - bool
            flag to set the Scale parameter free, NOTE: it is not
            recommended to set this flag to True
        BreakValue_free - bool
            flag to set the BreakValue parameter free
        Index2_free - bool
            flag to set the Index2 parameter free
        Beta_free - bool
            flag to set the Beta parameter free
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if the Prefactor, BreakValue, or Scale parameters are invalid,
            raise an error
        '''
        

        #do some sanity checks on inputs
        if Prefactor<0 or BreakValue<=0 or Scale<=0:
            raise ValueError("Input parameters invalid, at least one of Prefactor, BreakValue, or Scale is negative or zero.")

        #massage the inputs a bit and make
        #them accessible later if desired
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=Prefactor/self.Prefactor_scale
        self.Scale=Scale
        self.Index1=abs(Index1)
        self.BreakValue=BreakValue
        self.Index2=abs(Index2)
        self.Beta=Beta

        #set some descriptive attributes
        self.model_name='SmoothBrokenPowerLaw'
        self.norm_par='Prefactor'

        #create a dictionary for each parameter
        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','minimum':0,
                     'maximum':1e6,'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(Index1_free),'name':'Index1','minimum':0,
                     'maximum':10,'scale':-1,'value':self.Index1})

        self.parameters.append({'free':int(Scale_free),'name':'Scale','minimum':min(100,self.Scale*0.9),
                           'maximum':max(500000,self.Scale*1.5),'scale':1,'value':self.Scale})

        self.parameters.append({'free':int(BreakValue_free),'name':'BreakValue',
                     'minimum':min(100,self.BreakValue*0.9),
                     'maximum':max(500000,self.BreakValue*1.5),'scale':1,'value':self.BreakValue})

        self.parameters.append({'free':int(Index2_free),'name':'Index2','minimum':0,
                     'maximum':10,'scale':-1,'value':self.Index2})

        self.parameters.append({'free':int(Beta_free),'name':'Beta','minimum':0.01,
                     'maximum':10,'scale':1,'value':self.Beta})

        #now build the spectrum document element
        self.build()

    def ExpCutoff(self,Prefactor=1e-12,Index=2,Scale=1000,Ebreak=1000,
                  P1=100,P2=0,P3=0,Prefactor_free=True,Index_free=True,
                  Ebreak_free=True,P1_free=False,P2_free=False,P3_free=False,
                  Scale_free=False,model=None):

        '''
        create parameter dictionaries for ExpCutoff model and then
        build spectrum document element

        Parameters
        ----------
        Prefactor - float
            value of Prefactor parameter, normalization in MeV**-1 cm**-2 s**-1
        Index - float
            value of Index parameter, low-energy spectral index using E**-(Index) convention
        Scale - float
            value of Scale parameter, scale energy in MeV
        Ebreak - float
            value of Ebreak parameter, spectral break energy in MeV
        P1 - float
            value of the P1 parameter
        P2 - float
            value of the P2 parameter
        P3 - float
            value of the P3 parameter
        Prefactor_free - bool
            flag to set the Prefactor parameter free
        Index_free - bool
            flag to set the Index parameter free
        Scale_free - bool
            flag to set the Scale parameter free, NOTE: it is not
            recommended to set this flag to True
        EbreakV_free - bool
            flag to set the BreakValue parameter free
        P1_free - bool
            flag to set the P1 parameter free
        P2_free - bool
            flag to set the P2 parameter free
        P3_free - bool
            flag to set the P3 parameter free
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if Prefactor, Scale, or Ebreak inputs are invalid, raise
            an error
        '''
        
        #some checks on the inputs
        if Prefactor<0 or Scale<=0 or Ebreak<=0:
            raise ValueError("Input parameters invalid, at least one of Prefactor, Scale, or Ebreak is negative or zero.")

        #massage the inputs a bit and then make them
        #accessible later, if desired
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=Prefactor/self.Prefactor_scale
        self.Index=abs(Index)
        self.Scale=Scale
        self.Ebreak=Ebreak
        self.P1=P1
        self.P2=P2
        self.P3=P3

        #set more attributes
        self.model_name='ExpCutoff'
        self.norm_par='Prefactor'

        #create a dictionary for each parameters
        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','minimum':0,
                     'maximum':1e6,'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(Index_free),'name':'Index','minimum':0,
                     'maximum':10,'scale':-1,'value':self.Index})

        self.parameters.append({'free':int(Scale_free),'name':'Scale',
                     'minimum':min(100,self.Scale*0.9),
                     'maximum':max(500000,self.Scale*1.5),'scale':1,'value':self.Scale})

        self.parameters.append({'free':int(Ebreak_free),'name':'Ebreak',
                     'minimum':min(100,self.Ebreak*0.9),
                     'maximum':max(500000,self.Ebreak*1.5),'scale':1,'value':self.Ebreak})

        self.parameters.append({'free':int(P1_free),'name':'P1','minimum':0.1,
                     'maximum':max(300,self.P1*1.5),'scale':1,'value':self.P1})

        self.parameters.append({'free':int(P2_free),'name':'P2','minimum':-1,'maximum':1,
                     'scale':1,value:self.P2})

        self.parameters.append({'free':int(P3_free),'name':'P3','minimum':-1,'maximum':1,
                     'scale':1,value:self.P3})

        #now build the spectrum document element
        self.build()

    def BPLExpCutoff(self,Prefactor=1e-12,Index1=1.8,Index2=2.3,BreakValue=1000,
                  Eabs=10,P1=100,Prefactor_free=True,Index1_free=True,Index2_free=False,
                  BreakValue_free=True,Eabs_free=True,P1_free=False,model=None):
        '''
        create parameter dictionaries for BPLExpCutoff model and then build spectrum
        document element

        Parameters
        ----------
        Prefactor - float
            value of the Prefactor parameter, units of MeV**-1 cm**-2 s**-1
        Index1 - float
            value of the Index1 parameter, low-energy spectral index using E**(-Index1) convention
        Index2 - float
            value of the Index2 parameter, high-energy spectral index using E**(-Index2) convention
        BreakValue - float
            value of the BreakValue parameter, break energy in units of MeV
        Eabs - float
            value of the Eabs parameter
        P1 - float
            value of the P1 parameters
        Prefactor_free - bool
            flag to set the Prefactor parameter free
        Index1_free - bool
            flag to set the Index1 parameter free
        Index2_free - bool
            flag to set the Index2 parameter free
        BreakValue_free - bool
            flag to set the BreakValue parameter free
        Eabs_free - bool
            flag to set the Eabs parameter free
        P1_free - bool
            flag to set the P1 parameter free
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if Prefactor, BreakValue, or Eabs inputs are invalid, raise an error
        '''
        
        #some checks on the inputs
        if Prefactor<0 or BreakValue<=0 or Eabs<=0:
            raise ValueError("Input parameters invalid, at least one of Prefactor, BreakValue, or Eabs is negative or zero.")

        #massage the inputs a bit and make them
        #accessible later if needed
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=Prefactor/self.Prefactor_scale
        self.Index1=abs(Index1)
        self.Index2=Index2
        self.BreakValue=BreakValue
        self.Eabs=Eabs
        self.P1=P1

        #set some more attributes
        self.model_name='BPLExpCutoff'
        self.norm_par='Prefactor'

        #make a dictionary for each parameter
        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','minimum':0,
                     'maximum':1e6,'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(Index1_free),'name':'Index1','minimum':0,
                     'maximum':10,'scale':-1,'value':self.Index1})

        self.parameters.append({'free':int(Index2_free),'name':'Index2','minimum':0,
                     'maximum':10,'scale':-1,'value':self.Index2})

        self.parameters.append({'free':int(BreakValue_free),'name':'BreakValue',
                     'minimum':min(100,self.BreakValue*0.9),
                     'maximum':max(500000,self.BreakValue*1.5),'scale':1,'value':self.BreakValue})

        self.parameters.append({'free':int(Eabs_free),'name':'Eabs',
                     'minimum':min(1,self.Eabs*0.9),
                     'maximum':max(300,self.Eabs*1.5),'scale':1,'value':self.Eabs})

        self.parameters.append({'free':int(P1_free),'name':'P1','minimum':0.1,
                     'maximum':max(300,self.P1*1.5),'scale':1,'value':self.P1})

        #build the spectrum document element
        self.build()

    def PLSuperExpCutoff(self,Prefactor=1e-12,Index1=2,Scale=1000,Cutoff=1000,Index2=1,
                         Prefactor_free=True,Index1_free=True,Cutoff_free=True,Index2_free=False,
                         Scale_free=False,model=None):
        '''
        create parameter dictionaries for the PLSuperExpCutoff model and then build
        spectrum parameter element

        Parameters
        ----------
        Prefactor - float
            value of the Prefactor parameter, in units of MeV**-1 cm**-2 s**-1
        Index1 - float
            value of the Index1 parameter, uses the E**(-Index1) convention
        Scale - float
            value of the Scale parameter, in units of MeV
        Cutoff - float
            value of the Cutoff parameter, in units of MeV
        Index2 - float
            value of the Index2 parameter (the exponential index)
        Prefactor_free - bool
            flag to set the Prefactor parameter free
        Index1_free - bool
            flag to set the Index1 parameter free
        Cutoff_free - bool
            flag to set the Cutoff parameter free
        Index2_free - bool
            flag to set the Index2 parameter free
        Scale_free - bool
            flag to set the Scale parameter free, NOTE: it is not
            recommended to set this flag to True
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if Prefactor, Scale, Cutoff, or Index2 inputs are invalid, raise
            an error
        '''
        
        #some checks on the inputs
        if Prefactor<0 or Scale<=0 or Cutoff<=0 or Index2<0:
            raise ValueError("Input parameters invalid, at least one of Prefactor, Scale, Cutoff, or Index2 is negative or zero.")

        #massage the inputs a little and make
        #them available later if desired
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=Prefactor/self.Prefactor_scale
        self.Index1=abs(Index1)
        self.Scale=Scale
        self.Cutoff=Cutoff
        self.Index2=Index2

        #set some more attributes
        self.model_name='PLSuperExpCutoff'
        self.norm_par='Prefactor'

        #make a dictionary for each parameter
        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','minimum':0,
                     'maximum':1e6,'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(Index1_free),'name':'Index1','minimum':0,
                     'maximum':10,'scale':-1,'value':self.Index1})

        self.parameters.append({'free':int(Scale_free),'name':'Scale','minimum':min(100,self.Scale*0.9),
                     'maximum':max(500000,self.Scale*1.5),'scale':1,'value':self.Scale})

        self.parameters.append({'free':int(Cutoff_free),'name':'Cutoff','minimum':min(100,self.Cutoff*0.9),
                     'maximum':max(500000,self.Cutoff*1.5),'scale':1,'value':self.Cutoff})

        self.parameters.append({'free':int(Index2_free),'name':'Index2','minimum':0,'maximum':max(5,self.Index2*1.5),
                     'scale':1,'value':self.Index2})

        #create the spectrum document element
        self.build()

    def PLSuperExpCutoff2(self,Prefactor=1e-12,Index1=2,Scale=1000,Expfactor=1,Index2=1,
                          Prefactor_free=True,Index1_free=True,Expfactor_free=True,Index2_free=False,
                          Scale_free=False,model=None):
        '''
        create parameter dictionaries for the PLSuperExpCutoff2 model and then build
        spectrum parameter element

        Parameters
        ----------
        Prefactor - float
            value of the Prefactor parameter, in units of MeV**-1 cm**-2 s**-1
        Index1 - float
            value of the Index1 parameter, uses the E**(-Index1) convention
        Scale - float
            value of the Scale parameter, in units of MeV
        Expfactor - float
            value of the Expfactor parameter
        Index2 - float
            value of the Index2 parameter (the exponential index)
        Prefactor_free - bool
            flag to set the Prefactor parameter free
        Index1_free - bool
            flag to set the Index1 parameter free
        Expfactor_free - bool
            flag to set the Cutoff parameter free
        Index2_free - bool
            flag to set the Index2 parameter free
        Scale_free - bool
            flag to set the Scale parameter free, NOTE: it is not
            recommended to set this flag to True
        model - None-type or str
            model name, not used, included only for **kwargs compatability

        Raises
        ------
        ValueError
            if Prefactor, Scale, or Index2 inputs are invalid, raise
            an error
        '''
        
        #some checks on the inputs
        if Prefactor<0 or Scale<=0 or Index2<0:
            raise ValueError("Input parameters invalid, at least one of Prefactor, Scale, or Index2 is negative or zero.")

        #massage the inputs a little and assign
        #to attributes for access later if desired
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=Prefactor/self.Prefactor_scale
        self.Index1=abs(Index1)
        self.Scale=Scale
        self.Expfactor=Expfactor/1e-3
        self.Index2=Index2

        #set a couple more, useful attributes
        self.model_name='PLSuperExpCutoff2'
        self.norm_par='Prefactor'

        #make a dictionary for each parameter
        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','minimum':0,
                     'maximum':1e6,'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(Index1_free),'name':'Index1','minimum':0,
                     'maximum':10,'scale':-1,'value':self.Index1})

        self.parameters.append({'free':int(Scale_free),'name':'Scale','minimum':min(100,self.Scale*0.9),
                     'maximum':max(500000,self.Scale*1.5),'scale':1,'value':self.Scale})

        self.parameters.append({'free':int(Expfactor_free),'name':'Expfactor','minimum':-1,
                     'maximum':max(100,self.Expfactor*1.5),'scale':1e-3,'value':self.Expfactor})

        self.parameters.append({'free':int(Index2_free),'name':'Index2','minimum':0,'maximum':max(5,self.Index2*1.5),
                     'scale':1,'value':self.Index2})

        #create the spectrum document element
        self.build()

    def PLSuperExpCutoff3(self,Prefactor=1e-12,IndexS=2,Scale=1000,Expfactor2=2,Index2=1,
                          Prefactor_free=True,IndexS_free=True,Expfactor2_free=True,Index2_free=False,
                          Scale_free=False,model=None):
        '''
        create parameter dictionaries for the PLSuperExpCutoff3 model and then build
        spectrum parameter element

        Parameters
        ----------
        Prefactor - float
            value of the Prefactor parameter, in units of MeV**-1 cm**-2 s**-1
        IndexS - float
            value of the IndexS parameter, uses the E**(-IndexS) convention
        Scale - float
            value of the Scale parameter, in units of MeV
        Expfactor2 - float
            value of the Cutoff parameter, in units of MeV
        Index2 - float
            value of the Index2 parameter (the exponential index)
        Prefactor_free - bool
            flag to set the Prefactor parameter free
        IndexS_free - bool
            flag to set the Index1 parameter free
        Expfactor2_free - bool
            flag to set the Cutoff parameter free
        Index2_free - bool
            flag to set the Index2 parameter free
        Scale_free - bool
            flag to set the Scale parameter free, NOTE: it is not
            recommended to set this flag to True
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if Prefactor, Scale, or Index2 inputs are invalid, raise
            an error
        '''
        
        #some checks on the inputs
        if Prefactor<0 or Scale<=0 or Index2<0:
            raise ValueError("Input parameters invalid, at least one of Prefactor, Scale, or Index2 is negative or zero.")

        #massage the inputs a bit and assign to
        #attributes for access later if desired
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=Prefactor/self.Prefactor_scale
        self.IndexS=abs(IndexS)
        self.Scale=Scale
        self.Expfactor2=Expfactor2
        self.Index2=Index2

        #make some additional, useful attributes
        self.model_name='PLSuperExpCutoff3'
        self.norm_par='Prefactor'

        #create a dictionary for each parameter
        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','minimum':0,
                     'maximum':1e6,'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(IndexS_free),'name':'IndexS','minimum':0,
                     'maximum':10,'scale':-1,'value':self.IndexS})

        self.parameters.append({'free':int(Scale_free),'name':'Scale','minimum':min(100,self.Scale*0.9),
                     'maximum':max(500000,self.Scale*1.5),'scale':1,'value':self.Scale})

        self.parameters.append({'free':int(Expfactor2_free),'name':'Expfactor2','minimum':0,
                     'maximum':10,'scale':1,'value':self.Expfactor2})

        self.parameters.append({'free':int(Index2_free),'name':'Index2','minimum':0,'maximum':max(5,self.Index2*1.5),
                     'scale':1,'value':self.Index2})

        #create the spectrum document element
        self.build()

    def PLSuperExpCutoff4(self,Prefactor=1e-12,IndexS=2,Scale=1000,ExpfactorS=1,Index2=1,
                          Prefactor_free=True,IndexS_free=True,ExpfactorS_free=True,Index2_free=False,
                          Scale_free=False,model=None):
        '''
        create parameter dictionaries for the PLSuperExpCutoff4 model and then build
        spectrum parameter element

        Parameters
        ----------
        Prefactor - float
            value of the Prefactor parameter, in units of MeV**-1 cm**-2 s**-1
        Indexs - float
            value of the IndexS parameter, uses the E**(-IndexS) convention
        Scale - float
            value of the Scale parameter, in units of MeV
        ExpfactorS - float
            value of the Cutoff parameter, in units of MeV
        Index2 - float
            value of the Index2 parameter (the exponential index)
        Prefactor_free - bool
            flag to set the Prefactor parameter free
        IndexS_free - bool
            flag to set the IndexS parameter free
        ExpfactorS_free - bool
            flag to set the Cutoff parameter free
        Index2_free - bool
            flag to set the Index2 parameter free
        Scale_free - bool
            flag to set the Scale parameter free, NOTE: it is not
            recommended to set this flag to True
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if Prefactor, Scale, or Index2 inputs is invalid, raise an error
        '''
        
        #some checks on the inputs
        if Prefactor<0 or Scale<=0 or Index2<0:
            raise ValueError("Input parameters invalid, at least one of Prefactor, Scale, or Index2 is negative or zero.")

        #massage the inputs a little and then assign to
        #attributes for later access if desired
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=Prefactor/self.Prefactor_scale
        self.IndexS=abs(IndexS)
        self.Scale=Scale
        self.ExpfactorS=ExpfactorS/1e-1
        self.Index2=Index2
        
        #assign a couple more useful attributes
        self.model_name='PLSuperExpCutoff4'
        self.norm_par='Prefactor'

        #create a dictionary for each parameter
        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','minimum':0,
                     'maximum':1e6,'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(IndexS_free),'name':'IndexS','minimum':0,
                     'maximum':10,'scale':-1,'value':self.IndexS})

        self.parameters.append({'free':int(Scale_free),'name':'Scale','minimum':min(100,self.Scale*0.9),
                     'maximum':max(500000,self.Scale*1.5),'scale':1,'value':self.Scale})

        self.parameters.append({'free':int(ExpfactorS_free),'name':'ExpfactorS','minimum':-10,
                     'maximum':max(100,self.ExpfactorS*1.5),'scale':1e-1,'value':self.ExpfactorS})

        self.parameters.append({'free':int(Index2_free),'name':'Index2','minimum':0,'maximum':max(5,self.Index2*1.5),
                     'scale':1,'value':self.Index2})

        #create the spectrum document element
        self.build()

    def LogParabola(self,norm=1e-9,alpha=1,beta=2,Eb=300,
                    norm_free=True,alpha_free=True,beta_free=True,Eb_free=False,model=None):
        '''
        create parameter dictionaries for the LogParabola model and then build
        spectrum document element

        Parameters
        ----------
        norm - float
            value of the norm paramter, in units of MeV**-1 cm**-2 s**-1
        alpha - float
            value of the alpha parameter
        beta - float
            value of the beta parameter
        Eb - float
            value of the Eb parameter, in units of MeV
        norm_free - bool
            flag to set the norm parameter free
        alpha_free - bool
            flag to set the alpha parameter free
        beta_free - bool
            flag to set the beta paramter free
        Eb_free - bool
            flag to set the Eb parameter free
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if norm, alpha, or Eb inputs are invalid, raise an error
        '''
        
        #some checks on the inputs
        if norm<0 or alpha<0 or Eb<=0:
            raise ValueError("Input parameters invalid, at least one of norm, alpha, or Eb is negative or zero.")

        #massage in the inputs a bit and then assign to
        #attributes for later access if desired
        self.norm_scale=10**np.floor(np.log10(norm))
        self.norm=norm/self.norm_scale
        self.alpha=alpha
        self.beta=beta
        self.Eb=Eb

        #asign some more useful attributes
        self.model_name='LogParabola'
        self.norm_par='norm'

        #create a dictionary for each parameter
        self.parameters=[{'free':int(norm_free),'name':'norm','minimum':0,
                     'maximum':1e6,'scale':self.norm_scale,'value':self.norm}]

        self.parameters.append({'free':int(alpha_free),'name':'alpha','minimum':0,
                     'maximum':max(5,self.alpha*1.5),'scale':1,'value':self.alpha})

        self.parameters.append({'free':int(beta_free),'name':'beta','minimum':min(-5,self.beta*1.5),
                     'maximum':max(10,self.beta*1.5),'scale':1,'value':self.beta})

        self.parameters.append({'free':int(Eb_free),'name':'Eb','minimum':min(10,self.Eb*0.9),
                           'maximum':max(500000,self.Eb*1.5),'scale':1,'value':self.Eb})

        #create the spectrum document element
        self.build()

    def Gaussian(self,Prefactor=1e-9,Mean=7e4,Sigma=1e3,Prefactor_free=True,
                 Mean_free=True,Sigma_free=True,model=None):
        '''
        create parameter dictionaries for the Gaussian model and then build
        spectrum document element

        Parameters
        ----------
        Prefactor - float
            value of the Prefactor parameter, units of cm**-2 s**-1
        Mean - float
            value of the Mean parameter, units of MeV
        Sigma - float
            value of the Sigma parameter, units of MeV
        Prefactor_free - bool
            flag to set the Prefactor parameter free
        Mean_free - bool
            flag to set the Mean parameter free
        Sigma_free - bool
            flag to set the Sigma parameter Free
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if Prefactor, Mean, or Sigma inputs are invalid, raise
            an error
        '''

        #some checks on the inputs
        if Prefactor<=0 or Mean<0 or Sigma<=0:
            raise ValueError("Input parameters invalid, at least one of Prefactor, Mean, or Sigma <=0")

        #massage the inputs a bit and then assign to
        #attributes for access later if desired
        self.Prefactor_scale=10**np.floor(np.log10(Prefactor))
        self.Prefactor=Prefactor/self.Prefactor_scale
        self.Mean=Mean
        self.Sigma=Sigma

        #assign a couple more, useful attributes
        self.model_name='Gaussian'
        self.norm_par='Prefactor'

        #create a dictionary for each parameter
        self.parameters=[{'free':int(Prefactor_free),'name':'Prefactor','minimum':0,'maximum':1e6,
                          'scale':self.Prefactor_scale,'value':self.Prefactor}]

        self.parameters.append({'free':int(Mean_free),'name':'Mean','minimum':min(1e2,self.Mean*0.9),
                          'maximum':max(1e5,self.Mean*1.5),'scale':1,'value':self.Mean})

        self.parameters.append({'free':int(Sigma_free),'name':'Sigma','minimum':min(30,self.Sigma*0.9),
                          'maximum':max(1e4,self.Sigma*1.5),'scale':1,'value':self.Sigma})

        #create the spectrum document element
        self.build()

    def ConstantValue(self,Value=1,Value_free=True,model=None):
        '''
        create parameter dictionaries for the ConstantValue model and then build the
        spectrum document element

        Parameters
        ----------
        Value - float
            value of the Value parameter
        Value_free - bool
            flag to set the Value parameter free
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if Value input is invalid, raise an error
        '''

        #quick check on the Value parameter
        if Value<0:
            raise ValueError(f"Invalue input for Value = {Value}, must be >=0.")

        #assign values to some attributes to be
        #available later if desired
        self.Value=Value

        self.model_name='ConstantValue'
        self.norm_par="Value"

        #create a dictionary for the parameter
        self.parameters=[{'free':int(Value_free),'name':'Value','minimum':0,
                          'maximum':max(10,self.Value*1.5),'scale':1,'value':self.Value}]

        #create the spectrum document element
        self.build()

    def BandFunction(self,norm=1e-9,alpha=1.8,beta=2.5,Ep=0.1,norm_free=True,
                     alpha_free=True,beta_free=True,Ep_free=True,model=None):
        '''
        create parameter dictionaries for BandFunction model and then build the
        spectrum document element

        Parameters
        ----------
        norm - float
            value of the norm parameter, untis of MeV**-1 cm**-2 s**-1
        alpha - float
            value of the alpha parameter
        beta - float
            value of the beta parameter
        Ep - float
            value of the Ep parameter
        norm_free - bool
            flag to set the norm parameter free
        alpha_free - bool
            flag to set the alpha parameter free
        beta_free - bool
            flag to set the beta parameter free
        Ep_free - bool
            flag to set the Ep parameter free
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if norm or Ep inputs are invalid, raise an error
        '''

        #some checks on the inputs
        if norm<0 or Ep<=0:
            raise ValueError('Invalid input parameters, at least one of norm or Ep is <=0')

        #massage the inputs a little and then assign to
        #attributes for access later if desired
        self.norm_scale=10**np.floor(np.log10(norm))
        self.norm=norm/self.norm_scale
        self.alpha=alpha
        self.beta=beata
        self.Ep=Ep

        #assign a couple more useful attributes
        self.model_name='BandFunction'
        self.norm_par='norm'

        #create a dictionary for each parameter
        self.parameters=[{'free':int(norm_free),'name':'norm','minimum':0,'maximum':1e6,
                          'scale':self.norm_scale,'value':self.norm}]

        self.parameters.append({'free':int(alpha_free),'name':'alpha','minimum':-5,'maximum':-1,
                                'scale':1,'value':self.alpha})

        self.parameters.append({'free':int(beta_free),'name':'beta','minimum':-5,'maximum':-1,
                                'scale':1,'value':self.beta})

        self.parameters.append({'free':int(Ep_free),'name':'Ep','minimum':min(10,self.Ep*0.9),
                                'maximum':max(10000,self.Ep*1.5),'scale':1,'value':self.Ep})

        #creat the spectrum document element
        self.build()

    def DMFitFunction(self,spectrum_file=os.path.join('$(FERMI_DIR)','refdata','fermi',
                    'Likelihood','gammamc_dif.dat'),norm=5e20,sigmav=3e-26,mass=10,
                      bratio=1,channel0=4,channel1=1,norm_free=True,sigmav_free=True,
                      mass_free=True,bratio_free=False,channel0_free=False,channel1_free=False,
                      model=None):
        '''
        create the parameter dictionaries for the DMFitFunction model and then build
        spectrum document element

        Parameters
        ----------
        spectrum_file - str
            full path to the DM spectrum file, in this case a string is desired
            over a Path-like
        norm - float
            value of the norm parameter
        sigmav - float
            value of the sigmav parameter
        mass - float
            value of the mass parameter
        bratio - float
            value of the bratio parameter
        channel0 - float
            value of the channel0 parameter
        channel1 - float
            value of the channel1 parameter
        norm_free - bool
            flag to set the norm parameter free
        sigmav_free - bool
            flag to set the sigmav parameter free
        mass_free - bool
            flag to set the mass parameter free
        bratio_free - bool
            flag to set the bratio parameter free
        channel0_free - bool
            flag to set the channel0 parameter free
        channel1_free - bool
            flag to set the channel1 parameter free
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        ValueError
            if norm input is invalid, raise an error
        '''

        self.spectrum_file=str(spectrum_file)

        #quick check on inupts
        if norm<0:
            raise ValueError("Invalid input value of {norm = }, must be >=0")

        #massage the inputs a bit and then assign to some
        #attributes for access later if desired
        self.norm_scale=10**np.floor(np.log10(norm))
        self.norm=norm/self.norm_scale
        self.sigmav_scale=10**np.floor(np.log10(sigmav))
        self.sigmav=sigmav/self.sigmav_scale
        self.mass=mass
        self.bratio=bratio
        self.channel0=channel0
        self.channel1=channel1

        #assign some more useful attributes
        self.model_name='DMFitFunction'
        self.norm_par='norm'

        #create a dictionary for each parameter
        self.parameters=[{'free':int(norm_free),'name':'norm','minimum':1e-5,'maximum':1e-5,
                          'scale':self.norm_scale,'value':self.norm}]

        self.parameters.append({'free':int(sigmav_free),'name':'sigmav','minimum':0,'maximum':5000,
                                'scale':self.sigmav_scale,'value':self.sigmav})

        self.parameters.append({'free':int(mass_free),'name':'mass','minimum':1,'maximum':5000,
                                'scale':1,'value':self.mass})

        self.parameters.append({'free':int(bratio_free),'name':'bratio','minimum':0,'maximum':1,
                                'scale':1,'value':self.bratio})

        self.parameters.append({'free':int(channel0_free),'name':'channel0','minimum':min(1,self.channel0*0.9),
                                'maximum':max(10,self.channel0*1.5),'scale':1,'value':self.channel0})

        self.parameters.append({'free':int(channel1_free),'name':'channel1','minimum':min(1,self.channel1*0.9),
                                'maximum':max(10,self.channel1*1.5),'scale':1,'value':self.channel1})

        #create the spectrum document element
        self.build()

        #make sure to set the extra 'file' attribute of the document element
        self.spectrum.setAttribute('file',self.spectrum_file)
        

    def FileFunction(self,spectrum_file,Normalization=1,apply_edisp='false',
                     Normalization_free=True,model=None):
        '''
        create parameter dictionary for the FileFunction model and then build
        the spectrum document element

        Parameters
        ----------
        spectrum_file - str
            full path to the spectrum file, in this case a string is desired over
            a Path-like
        Normalization - float
            value of the normalization parameter
        apply_edisp - str or bool or int or float
            value of apply_edisp attribute, ideally this should be "true" or "false",
            case insensitive, but logic is implemented for bool, int, and float as well,
            if input is of the latter two types, >0 means 'true' and <=0 means 'false'
        Normalization_free - bool
            flag to set the Normalization parameter free
        model - None-type or str
            model name, not used, included only for **kwargs compatability
        
        Raises
        ------
        TypeError
            if apply_edisp input is not str, bool, int, or float type, raise an error
        ValueError
            if Normalization input is invalid, raise an error
            if apply_edisp is a str but not 'true' or 'false' (case insensitive),
              raise an error
        '''
        
        #some checks on the inputs
        if Normalization<0:
            raise ValueError(f"Input value of {Normalization = } is invalid, must be >0.")

        #try to allow for flexibility of the apply_edisp parameter
        if isinstance(apply_edisp,str):
            if apply_edisp.lower() not in ['true','false']:
                raise ValueError(f'Invalid input of {apply_edisp = }, must be "true" or "false".')

        elif isinstance(apply_edisp,bool):
            apply_edisp=str(apply_edisp).lower()

        elif isinstance(apply_edisp,(int,float)):
            apply_edisp='true' if apply_edisp>0 else 'false'

        else:
            raise TypeError(f'Invalid input type of {type(apply_edisp)}, must be str, bool, int, or float.')

        #assign values to some potentially helpful attributes
        self.Normalization=Normalization
        self.spectrum_file=str(spectrum_file)
        self.apply_edisp=apply_edisp.lower()

        self.model_name='FileFunction'
        self.norm_par='Normalization'

        #create a parameter dictionary
        self.parameters=[{'free':int(Normalization_free),'name':'Normalization',
                    'minimum':min(0.01,self.Normalization*0.9),'maximum':max(10,self.Normalization*1.5),
                    'scale':1,'value':self.Normalization}]

        #create the spectrum document element
        self.build()

        #add the extra attributes of the document element
        self.spectrum.setAttribute('apply_edisp',self.apply_edisp)
        self.spectrum.setAttribute('file',self.spectrum_file)

###########################################################################################
#class for the spatial models
#will only need to add the spatial class attribute/object to a larger XML document
###########################################################################################

class Spatial:
    '''
    A class to create a spectrum XML document element to be added to a source
    document element in a Fermi LAT XML region model

    ...

    Attributes
    ----------
    document : XML document
        document instance containing the spectrum document element
    functions : dict
        dictionary with keys as spectral function names and
        values pointing to the corresponding functions)
    spatial : XML document element
        document element corresponding to the spatialModel
    spatial_model : str
        name of the spatial model

    Methods
    -------
    get_function_dictionary()
        create the functions attribute as a dictionary of spectral model
        functions
    
    <for descriptions of all spatial models referenced in the methods below
    see https://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/source_models.html>
    ConstantValue(Value=1,spatial_model=None)
        create spatialModel document element corresponding to the ConstantValue
        spatial model
    MapCubeFunction(spatial_file,Normalization=1,spatial_model=None)
        create spatialModel document element corresponding to the MapCubeFunction
        spatial model
    RadialDisk(RA,DEC,Radius,spatial_model=None)
        create spatialModel document element corresponding to the RadialDisk
        spatial model
    RadialGaussian(RA,DEC,Sigma,spatial_model=None)
        create spatialModel document element corresponding to the RadialGaussian
        spatial model
    SkyDir(RA,DEC,spatial_model=None)
        create spatialModel document element corresponding to the SkyDir
        spatial model
    SpatialMap(spatial_file,Prefactor=1,spatial_model=None)
        create spatialModel document element corresponding to the SpatialMap
        spatial model
    '''
    
    def __init__(self,spatial_model,**kwargs):
        '''
        initialize the object, creating a XML document and document element
        then build the corresponding spatialModel element

        Parameters
        ----------
        spatial_model - str
            name of the spatial model
        **kwargs - dict or other keyword arguments
            keyword arguments corresponding to the necessary parameters of
            the input spatial_model
        
        Raises
        ------
        ValueError
            if requested spatial_model is not in the function dictionary,
            raise an error
        '''

        #create document and document element
        self.document=minidom.getDOMImplementation().createDocument(None,None,None)
        self.spatial=self.document.createElement('spatialModel')

        #assign the functions attribute
        self.get_function_dictionary()

        #do some checks on the input
        if spatial_model not in self.functions.keys():
            raise ValueError('{spatial_model = } not implemented')
        
        self.spatial_model=spatial_model
        
        #based on the input spatial_model
        self.functions[spatial_model](**kwargs)

    def get_function_dictionary(self):
        '''
        construct functions attribute as a dictionary
        of available spatial model functions
        '''
        
        self.functions={'SkyDir':self.SkyDir,
                   'RadialDisk':self.RadialDisk,
                   'RadialGaussian':self.RadialGaussian,
                   'SpatialMap':self.SpatialMap,
                   'ConstantValue':self.ConstantValue,
                   'MapCubeFunction':self.MapCubeFunction}

    def SkyDir(self,RA,DEC,spatial_model=None):
        '''
        assign necessary attributes and create document element
        children corresponding to the SkyDir spatial model

        Paramters
        ---------
        RA - float or str
            right ascension (J2000) of the source, in degrees
        DEC - float or str
            declination (J2000) of the source, in degrees
        spatial_model - None-type or str
            name of the spatial model, not used but necessary for
            kwargs functionality
        
        Raises
        ------
        ValueError
            if input RA or DEC value is outside of valid ranges, raise
            an error
        '''
        
        #do some sanity checks on the inputs
        if abs(RA)>360:
            raise ValueError(f'Input of {RA = } is invalid, must be between -360 and +360.')
        if abs(DEC)>90:
            raise ValueError(f'Input of {DEC = } is invalid, must be between -90 and 90.')

        #assign attributes for possible use later
        self.RA=RA
        self.DEC=DEC

        #set the type attribute
        self.spatial.setAttribute('type','SkyDirFunction')

        #create the document children
        self.spatial.appendChild(parameter_element(free='0',name='RA',maximum='360.0',
                minimum='-360.0',scale='1.0',value=f'{self.RA:.4f}'))
        self.spatial.appendChild(parameter_element(free='0',name='DEC',maximum='90.0',
                minimum='-90.0',scale='1.0',value=f'{self.DEC:.4f}'))

    
    def RadialDisk(self,RA,DEC,Radius,spatial_model=None):
        '''
        assign necessary attributes and create document element
        children corresponding to the RadialDisk spatial model

        Parameters
        ----------
        RA - float or str
            right ascension (J2000) of source in degrees
        DEC - float or str
            declination (J2000) of source in degrees
        Radius - float or str
            radius of extended disk in degrees
        spatial_model - None-type or str
            name of the spatial model, not used but necessary for
            kwargs functionality
        
        Raises
        ------
        ValueError
            if RA, DEC, or Radius inputs are invalid, raise an error
        '''

        #do some sanity checks, including converting Radius to float
        if abs(float(RA))>360:
            raise ValueError(f'Input value {RA = } is invalid, must be between -360 and +360.')

        if abs(float(DEC))>90:
            raise ValueError(f'Input value {DEC = } is invalid, must be between -90 and 90.')

        if float(Radius)<=0:
            raise ValueError(f'Input value {Radius = } is invalid, must be >0.')

        #save inputs to attributes for possible
        #use later
        self.RA=float(RA)
        self.DEC=float(DEC)
        self.Radius=float(Radius)

        #set the type attribute
        self.spatial.setAttribute('type','RadialDisk')

        #create the necessary document children
        self.spatial.appendChild(parameter_element(free='0',name='Radius',
                maximum=max(10,self.Radius*1.5),minimum='0',scale='1',value=self.Radius))
        self.spatial.appendChild(parameter_element(free='0',name='RA',maximum='360.0',
                minimum='-360.0',scale='1.0',value=f'{self.RA:.4f}'))
        self.spatial.appendChild(parameter_element(free='0',name='DEC',maximum='90.0',
                minimum='-90.0',scale='1.0',value=f'{self.DEC:.4f}'))

    def RadialGaussian(self,RA,DEC,Sigma,spatial_model=None):
        '''
        assign necessary attributes and create document element
        children corresponding to the RadialGaussian sptial model

        Paramters
        ---------
        RA - float or str
            right ascension (J2000) of source in degrees
        DEC - float or str
            declination (J2000) of source in degrees
        Sigma - float or str
            value of Sigma parmeter in degrees
        spatial_model - None-type or str
            name of the spatial model, not used but necessary for
            kwargs functionality
        
        Raises
        ------
        ValueError
            if RA, DEC, Sigma inputs are invalid, raise an error
        '''
        
        #do some sanity checks, including converting Sigma to a float
        if abs(float(RA))>360:
            raise ValueError(f'Input value of {RA = } is invalid, must be between -360 and +360.')
        if abs(float(DEC))>90:
            raise ValueError(f'Input value of {DEC = } is invalid, must be between -90 and 90.')
        if float(Sigma)<=0:
            raise ValueError(f'Input value of {Sigma = } is invalid, must be >0.')

        #assign inputs to attributes for possible
        #use later if desired
        self.RA=float(RA)
        self.DEC=float(DEC)
        self.Sigma=float(Sigma)

        #set the type attribute
        self.spatial.setAttribute('type','RadialGaussian')

        #create the document element children
        self.spatial.appendChild(parameter_element(free='0',name='Sigma',
                maximum=max(10,self.Sigma*1.5),minimum='0',scale='1',value=self.Sigma))
        self.spatial.appendChild(parameter_element(free='0',name='RA',maximum='360.0',
                minimum='-360.0',scale='1.0',value=f'{self.RA:.4f}'))
        self.spatial.appendChild(parameter_element(free='0',name='DEC',maximum='90.0',
                minimum='-90.0',scale='1.0',value=f'{self.DEC:.4f}'))

    def SpatialMap(self,spatial_file,Prefactor=1,spatial_model=None):
        '''
        assign necessary attributes and create document element
        children corresponding to the SpatialMap spatial model

        Parameters
        ----------
        spatial_file - str
            full path to the spatial_file, in this case a string is desired
            over a Path-like
        Prefactor - str or float
            value of Prefactor parameter, units of MeV**-1 sr**-1 cm**-2 s**-1
        spatial_model - None-type or str
            name of the spatial model, not used but necessary for
            kwargs functionality
        
        Raises
        ------
        ValueError
            if Prefactor input is invalid, raise an error
        '''

        #quick check on input
        if float(Prefactor)<=0:
            raise ValueError(f'Input of {Prefactor = :e} is invalid, must be >=0')

        #assign inputs to some attributes for access later
        self.spatial_file=str(spatial_file)
        self.Prefactor=Prefactor

        #assign the necessary attributes to the spatialModel document element
        self.spatial.setAttribute('type','SpatialMap')
        self.spatial.setAttribute('file',self.spatial_file)
        self.spatial.setAttribute('map_based_integral','true')

        #create the document element child
        self.spatial.appendChild(parameter_element(free="0",name="Prefactor",maximum="1000",
                minimum="0.001",scale="1",value=self.Prefactor))

    def ConstantValue(self,Value=1,spatial_model=None):
        '''
        assign necessary attributes and create document element
        children corresponding to the ConstantValue spatial model

        Parameters
        ----------
        Value - float or string

        Raises
        ------
        ValueError
            if Value input is invalid, raise an error
        '''
        
        #check on inputs
        if float(Value)<0:
            raise ValueError(f'Input {Value = } is invalid, must be >=0.')

        #assign input to some attribute for access later
        self.Value=float(Value)

        #set type attribute
        self.spatial.setAttribute('type','ConstantValue')

        #create document element child
        self.spatial.appendChild(parameter_element(free="0",name="Value",maximum=max(10,self.Value*1.5),
                    minimum='0',scale='1',value=self.Value))

    def MapCubeFunction(self,spatial_file,Normalization=1,spatial_model=None):
        '''
        assign necessary attributes and create document element
        children corresponding to the MapCubeFunction sptial model

        Paramters
        ---------
        spatial_file - str
            full path to the spatial_file file, in this case a string is
            desired over a Path-like
        Normalization - float or str
            value of the Normalization parameter
        spatial_model - None-type or str
            name of the spatial model, not used but necessary for
            kwargs functionality
        
        Raises
        ------
        ValueError
            if Normalization input is invalid, raise an error
        '''
        
        #check on input value
        if (Normalization)<0:
            raise ValueError(f'Input value of {Normalization = } is invalid, must be >0.')

        #assign inputs to attributes for possible access
        #later if desired
        self.spatial_file=str(spatial_file)
        self.Normalization=float(Normalization)

        #set the necessary attributes
        self.spatial.setAttribute('type','MapCubeFunction')
        self.spatial.setAttribute('file',self.spatial_file)

        #create the document element child
        self.spatial.appendChild(parameter_element(free="0",name="Normalization",
                maximum=max(1000,self.Normalization*1.5),
                minimum=min(0.0,self.Normalization*0.9),scale=1,value=self.Normalization))
