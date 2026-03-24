
'''
Simplified biophysical model for investigating the effect of AIS position.

This is a model with a spherical isopotential soma, a large dendrite and a myelinated axon. 
The AIS contains a higher density of Na and K channels than the soma, the dendrite and the rest of the axon.
The axon beyond 50 um from the soma is myelinated.

Based on Goethals & Brette, Elife. 2020 Mar 30:9:e53432. doi: 10.7554/eLife.53432. 

'''
from brian2 import *

__all__ = ['model_na_kv1_myelin']

def model_na_kv1_myelin(params, demyelinated=False, resting_vm = -75.*mV, AIS_center = 5.*um, \
                gna_tot = 350.*nS, gk_tot = 150.*nS, morpho = None):
    '''
    
    params: a file that contains  the passive and channels parameters,
    
    neuron.I_CC is a current-clamp stimulation.
    neuron.V_VC is a voltage-clamp stimulation.
    neuron.VC_on is a voltage-clamp switch (0/1).
    neuron.I_VC is the voltage-clamp current.
    '''

    ### Morphology of the neuron
    morpho = Soma(15.2* um)
    dend_diam = 2.*um 
    dend_length = 1000.*um 
    axon_diam = 0.5 * um
    axon_length = 500. * um
    dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=500)
    axon = Cylinder(diameter=axon_diam, length=axon_length, n=500)
    morpho.dendrite = dendrite
    morpho.axon = axon

    ### Passive parameters
    EL = resting_vm 
    Cm = params.Cm 
    gL = params.gL
    Ri = params.Ri 
    
    ### Myelin parameters
    if demyelinated:
        Cm_myelin = Cm
        gL_myelin = gL
    else:
        Cm_myelin = params.Cm/4
        gL_myelin = gL/4
 
    ### Capacitance distribution 
    Cm_neuron = hstack([Cm/(uF / cm ** 2) * ones(morpho.dendrite.n + 1 + 50), Cm_myelin/(uF / cm ** 2) * ones(morpho.axon.n - 50)]) * uF / cm ** 2

    # Na channels parameters
    ENa = params.ENa 
    Gna = gna_tot 
    Va_soma = params.Va_soma
    Va = params.Va #of the AIS
    Ka = params.Ka 
    Taum_max = params.Taum_max 

    Vh_soma = params.Vh_soma 
    Vh = params.Vh #of the AIS
    Kh = params.Kh 
    Tauh_max = params.Tauh_max

    # K channels parameters
    EK = params.EK 
    Gk = gk_tot

    Vn = params.Vn 
    Kn = params.Kn 
    Taun_max = params.Taun_max 

    # Na channels conductances
    gna_soma = params.gna_soma 
    gna_dend = params.gna_dend 

    # K channels parameters
    gk_soma = params.gk_soma 
    gk_dend = params.gk_dend 

    # Equations
    eqs = '''
    Im = (gL*(EL-v) + gNa*m*h*(ENa-v) + gK*n**8*(EK-v)) : amp/meter**2
    INa = gNa*m*h*(ENa-v) : amp/meter**2
    IK = gK*n**8*(EK-v) : amp/meter**2

    dm/dt = alpham*(1-m) - betam*m : 1
    dh/dt = alphah*(1-h) - betah*h : 1
    dn/dt = alphan*(1-n) - betan*n : 1

    alpham = (1/ka)*(v-va) / (1-exp(-(v-va)/ka)) /(2*taum_max) : Hz
    betam = -(1/ka)*(v-va) / (1-exp((v-va)/ka)) /(2*taum_max) : Hz

    alphah = -(1/kh)*(v-vh) / (1-exp((v-vh)/kh)) /(2*tauh_max) : Hz
    betah = (1/kh)*(v-vh) / (1-exp(-(v-vh)/kh)) /(2*tauh_max) : Hz

    alphan = (1/kn)*(v-vn) / (1-exp(-(v-vn)/kn)) /(2*taun_max) : Hz
    betan = -(1/kn)*(v-vn) / (1-exp((v-vn)/kn)) /(2*taun_max): Hz
    
    gL: siemens/meter**2
    gNa : siemens/meter**2
    gK : siemens/meter**2
    va : volt
    vh : volt
    vn : volt
    ka : volt
    kh : volt
    kn : volt
    taum_max : second
    tauh_max : second
    taun_max : second

    I : amp (point current)
    '''
    
    # Current-clamp and voltage-clamp stimulation
    gclamp = 100*usiemens

    eqs += '''
    I_CC : amp (point current) 
    I_VC = gclamp*VC_on*(V_VC - v) : amp (point current)
    V_VC : volt
    VC_on : 1 # if 1, voltage-clamp is on
    
    I_hyp : amp (point current)
    '''

    neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=Cm_neuron, Ri=Ri, threshold="v > 0 *mV",  refractory=5*ms,
                                    namespace={'EL': EL, 'ENa': ENa, 'EK': EK,
                                                   'Ka': Ka, 'Va': Va, 'Taum_max': Taum_max,
                                                   'Kh': Kh, 'Vh': Vh, 'Tauh_max': Tauh_max,
                                                   'Kn': Kn, 'Vn': Vn, 'Taun_max': Taun_max,
                                                   'Va_soma': Va_soma, 
                                                   'Vh_soma': Vh_soma,                                                  
                                                   'gclamp': gclamp},
                                                   method="exponential_euler")

    # Parameters of the soma are assigned in the entire neuron for integration purpose, 
    # but gNa and gK at the soma and AIS are then modified
    neuron.va = Va_soma
    neuron.ka = Ka
    neuron.taum_max = Taum_max 

    neuron.vh = Vh_soma
    neuron.kh = Kh
    neuron.tauh_max = Tauh_max

    neuron.vn = Vn
    neuron.kn = Kn
    neuron.taun_max = Taun_max 

    #Dendrites and axon
    neuron.gNa = gna_dend
    neuron.gK = gk_dend
    neuron.gL = gL

    # Soma
    neuron.gNa[0] = gna_soma 
    neuron.gK[0] = gk_soma 
    neuron.gL[0] = gL
    
    # Initial segment
    initial_segment = morpho.axon[AIS_center]
    neuron.gNa[initial_segment] = Gna/neuron.area[initial_segment]
    neuron.gK[initial_segment] = Gk/neuron.area[initial_segment]
    neuron.va[initial_segment] = Va
    neuron.vh[initial_segment] = Vh
    
    ## Myelin
    myel_start= 50.*um
    myelinated_axon = morpho.axon[myel_start:]
    neuron.gL[myelinated_axon] = gL_myelin
    
    # Initialisation
    neuron.v = EL
    neuron.h = 1
    neuron.m = 0
    neuron.n = 1

    return neuron