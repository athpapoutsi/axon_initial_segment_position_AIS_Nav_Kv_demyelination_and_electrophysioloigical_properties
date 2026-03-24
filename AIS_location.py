from brian2 import *
from joblib import Parallel, delayed

import params_model_description
from model_analysis import (
    calculate_resting_state,
    firing_properties,
    measure_current_threshold,
    measure_input_resistance,
    measure_voltage_threshold,
)
from model_na_kv1_myelin import model_na_kv1_myelin

# Parameters
defaultclock.dt = 0.1*ms

params = params_model_description

GNa = 750.*nS # total Na conductance at the AIS
GK = 2800*nS # total K conductance at the AIS

centers = [27, 20]*um # AIS center positions
demyelinated=True

### SIMULATIONS 
# A function to measure electrophsyioogical properties with somatic current, for different AIS center positions
def BIO_model_AIS_position(ais_center, gna_tot, gk_tot):
    print ('AIS center:', ais_center, 'Gna:', gna_tot)
    params = params_model_description
    pulse_length = 10.*ms
    
    neuron = model_na_kv1_myelin(params=params, demyelinated=demyelinated, resting_vm =-71.*mV, AIS_center = ais_center, gna_tot=gna_tot, gk_tot=gk_tot)
    v_rest=calculate_resting_state(neuron)  

    neuron = model_na_kv1_myelin(params=params, demyelinated=demyelinated, resting_vm =-71.*mV, AIS_center = ais_center, gna_tot=gna_tot, gk_tot=gk_tot)
    i_rheo = measure_current_threshold(neuron, -71.*mV, pulse_length)

    neuron = model_na_kv1_myelin(params=params, demyelinated=demyelinated, resting_vm =-71.*mV, AIS_center = ais_center, gna_tot=gna_tot, gk_tot=gk_tot)
    v_thr= measure_voltage_threshold(neuron, ais_center, i_rheo, pulse_length)

    neuron = model_na_kv1_myelin(params=params, demyelinated=demyelinated, resting_vm =-71.*mV, AIS_center = ais_center, gna_tot=gna_tot, gk_tot=gk_tot)
    steady_pulse= measure_input_resistance(neuron, current=0.01*nA, resting_vm=-71.*mV)

    neuron = model_na_kv1_myelin(params=params, demyelinated=demyelinated, resting_vm =-71.*mV, AIS_center = ais_center, gna_tot=gna_tot, gk_tot=gk_tot)
    ff_prop= firing_properties(neuron, current=0.04*nA, resting_vm=-71.*mV)

    return v_rest, i_rheo, v_thr, steady_pulse, ff_prop

# run the simulations with joblib
if __name__ == '__main__':
    results = Parallel(n_jobs = 5)(delayed(BIO_model_AIS_position)(ais_center, GNa, GK) for ais_center in centers)
    
    v_rests=[results[k][0] for k in range(len(results))]
    i_thresholds= [results[k][1] for k in range(len(results))]
    v_thresholds= [results[k][2] for k in range(len(results))]

    IR=[results[k][3][0] for k in range(len(results))]
    voltage_IR= [results[k][3][1] for k in range(len(results))]
    time_IR= [results[k][3][2] for k in range(len(results))]

    ff=[results[k][4][0] for k in range(len(results))]
    voltage_ff= [results[k][4][1] for k in range(len(results))]
    time_ff= [results[k][4][2] for k in range(len(results))]
    v_ahps= [results[k][4][3] for k in range(len(results))]

print("\nAIS center, Vrest")
print(centers[0], v_rests[0])
print(centers[1], v_rests[1])
print("\nAIS center, I_rheobase, V_rheobase")
print(centers[0], i_thresholds[0], v_thresholds[0][0] )
print(centers[1], i_thresholds[1], v_thresholds[1][0])
print("\nAIS center, Input resistance")
print(centers[0], IR[0] )
print(centers[1], IR[1] )
print(" IR difference", IR[0]-IR[1])
print("\nAIS center, FF")
print(centers[0], ff[0]*2 )
print(centers[1], ff[1]*2 )
print("\nAIS center, AHP amplitude")
print(centers[0], v_thresholds[0][0]-v_ahps[0])
print(centers[1], v_thresholds[1][0]-v_ahps[1])


fig1=plt.figure()
plt.plot(time_IR[0], voltage_IR[0], 'k')
plt.plot(time_IR[1], voltage_IR[1], 'b')
plt.xlabel('time (ms)') 
plt.ylabel('Voltage (mV)') 
plt.title('Input resistance measurement')
# fig1.savefig(f"Myel_{GNa}_{GK}_IR.pdf", bbox_inches='tight', dpi=300)

fig2=plt.figure()
plt.plot(time_ff[0], voltage_ff[0], 'k')
plt.plot(time_ff[1], voltage_ff[1], 'b') 
plt.xlabel('time (ms)') 
plt.ylabel('Voltage (mV)') 
plt.title('Firing properties')
# # fig2.savefig(f"Myel_{GNa}_{GK}_ff.pdf", bbox_inches='tight', dpi=300)

fig3=plt.figure()
plt.plot(v_thresholds[0][2], v_thresholds[0][1], 'k')
plt.plot(v_thresholds[1][2], v_thresholds[1][1], 'b')
plt.plot(v_thresholds[0][2], v_thresholds[0][9], '0.5')
plt.plot(v_thresholds[1][2], v_thresholds[1][9], 'g')
plt.ylim([-91, -48])
plt.xlabel('time (ms)') 
plt.ylabel('Voltage (mV)') 
plt.title('Voltage threshold')
# # fig3.savefig(f"Myel_{GNa}_{GK}_VR.pdf", bbox_inches='tight', dpi=300)

fig4=plt.figure()
plt.plot(v_thresholds[0][2], v_thresholds[0][4], 'k')
plt.plot(v_thresholds[1][2], v_thresholds[1][4], 'b')
plt.plot(v_thresholds[0][2], v_thresholds[0][5], '0.5')
plt.plot(v_thresholds[1][2], v_thresholds[1][5], 'g')
plt.ylim([-1, 62])
plt.ylim([-0.1, 60])
plt.xlabel('time (ms)') 
plt.ylabel('Current (pA/μm^2)') 
plt.title('Sodium current at voltage threshold')
# fig4.savefig(f"Myel_{GNa}_{GK}_Ina.pdf", bbox_inches='tight', dpi=300)

fig5=plt.figure()
plt.plot([0, 10, centers[0]/um,40, 50], [max(v_thresholds[0][1]), max(v_thresholds[0][11]),max(v_thresholds[0][9]), max(v_thresholds[0][12]), max(v_thresholds[0][10])], 'k')
plt.plot([0, 10, centers[1]/um,40, 50], [max(v_thresholds[1][1]), max(v_thresholds[1][11]),max(v_thresholds[1][9]), max(v_thresholds[1][12]), max(v_thresholds[1][10])], 'b')
plt.xlabel('Distance from soma (um)')
plt.xticks([0, 10, centers[1]/um, centers[0]/um, 40, 50])
plt.ylabel('Voltage (mV)')
plt.ylim([-69, -48])
plt.title('Voltage gradient at voltage threshold')
# fig5.savefig(f"Myel_{GNa}_{GK}_V_distance.pdf", bbox_inches='tight', dpi=300)

plt.show()