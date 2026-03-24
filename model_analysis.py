
"""

Functions to measure electrical properties in neuron model.

"""
from brian2 import *

__all__ = ['calculate_resting_state', 'firing_properties', 'measure_current_threshold', 'measure_input_resistance', 'measure_voltage_threshold']

defaultclock.dt = 0.1*ms

def calculate_resting_state(neuron):
    '''
    Measure the resting membrane potential of a neuron model at the soma.
    '''
    latency=100*ms
    M = StateMonitor(neuron, ('v'), record = 0)
    run(latency+150.*ms)
        
    return M.v[0][int(latency/defaultclock.dt):].mean()/mV

def measure_input_resistance(neuron, current=0.01*nA, resting_vm=-75.*mV):
    '''
    Measure the input resistance at the soma, as the voltage gradient at the end of a current pulse divided by the amplitude of the current pulse.
    '''
    
    M = StateMonitor(neuron, ('v'), record = 0)
    
    neuron.V_VC[0] = resting_vm 
    neuron.VC_on[0] =0
    run(200*ms)
    neuron.VC_on[0] = 0
    neuron.I_CC[0] = current
    run(500.*ms)
    neuron.I_CC[0] = 0*amp 
    run(200*ms)
    
    Rin_soma = (M.v[0][int((200*ms+500*ms)/defaultclock.dt)]-M.v[0][int((199*ms)/defaultclock.dt)]) / current    
    
    return Rin_soma, M.v[0]/mV, M.t/ms
    
def measure_current_threshold(neuron, resting_vm=-72.*mV, pulse_length = 2.*ms, latency = 200.*ms, ):
    '''
    Measure the current threshold as the minimal somatic current that elicits a spike, for a specific duration of the current pulse,
    with the bissection method. 
    
    Returns the minimal somatic current amplitude that elicits a spike at the AIS.    
    '''
    i_max = 3.*nA
    i_min = 0.*nA
    i_current = 0.1*nA
    spike = False
    
    M = StateMonitor(neuron, ('v', 'm', 'I_VC'), record = 0)
        
    store()
    
    n_it = 0
    while True:
        restore()
        
        neuron.V_VC[0] = resting_vm 
        neuron.VC_on[0] = 0
        neuron.I_CC[0] = 0*amp 
        run(latency)
        neuron.VC_on[0] = 0
        neuron.I_CC[0] = i_current
        run(pulse_length)
        neuron.I_CC[0] = 0*amp 
        run(20*ms)
                 
        m_max_ais = max(M.m[0][int(latency/defaultclock.dt):])    
        if m_max_ais >= 0.5 and abs(i_current - i_min) <= 0.1*pA and spike == False :
            break
        if m_max_ais <= 0.5:
            i_min = i_current
            spike = False
        else: 
            i_max = i_current
            spike = True
        
        i_current = 0.5*i_max + 0.5*i_min
        
        n_it += 1
        
        if n_it > 25:
            print ('Too much iterations')
            break
                
    return i_current

def measure_voltage_threshold(neuron, AIS_center, i_rheo = 0.5*nA, pulse_length = 2.*ms, ):
    '''
    Measures the voltage threshold as the maximal membrane potential reached during the largest non-spiking situation.
    A current pulse of an amplitude just below the current threshold is injected at the soma.
    
    Returns: voltage threshold at the soma, 10 and 40 um away from the soma, at the AIS and at the start of the myelinated axon.

    '''
    print ('Measuring the voltage threshold with ireho', i_rheo*0.9985)
  
    M = StateMonitor(neuron, ('v','I_VC', 'INa'), record = 0)
    M_AIS = StateMonitor(neuron, ('v', 'INa'), record = neuron.morphology.axon[AIS_center])
    M_Myel = StateMonitor(neuron, ('v', 'INa'), record = neuron.morphology.axon[50.*um])
    M_10 = StateMonitor(neuron, ('v', 'INa'), record = neuron.morphology.axon[10.*um])
    M_40 = StateMonitor(neuron, ('v', 'INa'), record = neuron.morphology.axon[40.*um])

    neuron.VC_on[0] = 0
    neuron.I_CC[0] =0 *amp
    run(200.*ms)
    neuron.VC_on[0] = 0
    neuron.I_CC[0] = i_rheo * 0.9985
    run(pulse_length)
    neuron.I_CC[0] = 0*amp 
    run(20*ms)
    
    th_rheo_soma = max(M.v[0][int(180.*ms/defaultclock.dt):]) 
    v0_soma = mean(M.v[0][int(180.*ms/defaultclock.dt):int(200.*ms/defaultclock.dt)])

    return th_rheo_soma, M.v[0]/mV, M.t/ms, v0_soma, M.INa[0]/(amp/meter**2), M_AIS.INa[0]/(amp/meter**2), M_Myel.INa[0]/(amp/meter**2), M_10.INa[0]/(amp/meter**2), M_40.INa[0]/(amp/meter**2), M_AIS.v[0]/mV, M_Myel.v[0]/mV, M_10.v[0]/mV, M_40.v[0]/mV

def firing_properties(neuron, current=0.04*nA, resting_vm=-71.*mV):
    '''
    Measure the number of spikes for a given current pulse of 500 ms duration.
    '''
    
    M = StateMonitor(neuron, ('v','I_VC'), record = 0)
    spikes = SpikeMonitor(neuron, record = 0, variables='v')
    
    neuron.V_VC[0] = resting_vm 
    neuron.VC_on[0] =0
    run(200*ms)
    neuron.VC_on[0] = 0
    neuron.I_CC[0] = current
    run(500.*ms)
    neuron.I_CC[0] = 0*amp 
    run(200*ms)
    
    v_ahp = min(M.v[0][int((200*ms+40*ms)/defaultclock.dt):int((599*ms)/defaultclock.dt)])
    return spikes.count[0], M.v[0]/mV, M.t/ms, v_ahp


    
    
    
    
    
    
    
    
    
    
    
    

