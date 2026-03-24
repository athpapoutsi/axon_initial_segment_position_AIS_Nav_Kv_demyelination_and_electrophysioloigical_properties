from brian2 import *
from joblib import Parallel, delayed
from matplotlib.colors import ListedColormap

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

test_param='gks'
#If this is set to true, axon is completely demyelinated.
demyelinated=False

params = params_model_description
if test_param == 'gnas':
    gs = linspace (690, 770, 5)*nS # total Na conductances at the AIS
    GK = 2800*nS # total K conductance at the AIS
elif test_param == 'gks':
    gs = linspace (2800, 3600, 5)*nS # total K conductances at the AIS
    GNa = 750.*nS # total Na conductance at the AIS
else: 
    print('test_param should be gnas or gks')

print(gs)
centers = [27, 20]*um # AIS center positions
n1=len(gs)
n2=len(centers)

### SIMULATIONS 
# A function to measure electrophsyioogical properties with somatic current, for different AIS center positions and gs
def BIO_model_AIS_position(center, gna_tot, gk_tot):
    print ('AIS center:', center, 'Gna:', gna_tot, 'Gk:', gk_tot)

    params = params_model_description
    pulse_length = 10.*ms
    
    neuron = model_na_kv1_myelin(params=params, demyelinated=demyelinated, resting_vm =-71.*mV, AIS_center = center, gna_tot=gna_tot, gk_tot=gk_tot)
    v_rest=calculate_resting_state(neuron)  

    neuron = model_na_kv1_myelin(params=params, demyelinated=demyelinated, resting_vm =-71.*mV, AIS_center = center, gna_tot=gna_tot, gk_tot=gk_tot)
    i_rheo = measure_current_threshold(neuron, -71.*mV, pulse_length)

    neuron = model_na_kv1_myelin(params=params, demyelinated=demyelinated, resting_vm =-71.*mV, AIS_center = center, gna_tot=gna_tot, gk_tot=gk_tot)
    v_thr= measure_voltage_threshold(neuron, center, i_rheo, pulse_length)

    neuron = model_na_kv1_myelin(params=params, demyelinated=demyelinated, resting_vm =-71.*mV, AIS_center = center, gna_tot=gna_tot, gk_tot=gk_tot)
    steady_pulse= measure_input_resistance(neuron, current=0.01*nA, resting_vm=-71.*mV)

    neuron = model_na_kv1_myelin(params=params, demyelinated=demyelinated, resting_vm =-71.*mV, AIS_center = center, gna_tot=gna_tot, gk_tot=gk_tot)
    ff_prop= firing_properties(neuron, current=0.04*nA, resting_vm=-71.*mV)

    return v_rest, i_rheo, v_thr, steady_pulse, ff_prop

# run the simulations with joblib
if __name__ == '__main__':
    if test_param == 'gnas':
        results = Parallel(n_jobs = 5)(delayed(BIO_model_AIS_position)(center, gna_tot, GK) for center in centers for gna_tot in gs)
    elif test_param == 'gks':
        results = Parallel(n_jobs = 5)(delayed(BIO_model_AIS_position)(center, GNa, gk_tot) for center in centers for gk_tot in gs)
    else:
        print('test_param should be gnas or gks')   

    v_rests=[results[k][0] for k in range(len(results))]
    v_rests = np.array(v_rests).reshape((n2,n1))*1e3
    i_thresholds= [results[k][1] for k in range(len(results))]
    i_thresholds = np.array(i_thresholds).reshape((n2,n1))*1e12
    v_thresholds= [results[k][2][0] for k in range(len(results))]
    v_thresholds = np.array(v_thresholds).reshape((n2,n1))*1e3

    IR=[results[k][3][0] for k in range(len(results))]
    IR = np.array(IR).reshape((n2,n1))*1e-6
    ff=[results[k][4][0]*2 for k in range(len(results))]
    ff = np.array(ff).reshape((n2,n1))
    v_ahps= [results[k][4][3]-results[k][2][0] for k in range(len(results))]
    v_ahps = np.array(v_ahps).reshape((n2,n1))*1e3

### Plots
# Custom colormap
top = get_cmap('Oranges_r', 128)
bottom = get_cmap('Blues', 128)
newcolors = vstack((top(linspace(0.5, 1, 128)),
                       bottom(linspace(0.5, 1, 128))))
newcmp = ListedColormap(newcolors, name='OrangeBlue')

# Figure
cmap = get_cmap(newcmp)
colors = [cmap(i) for i in np.linspace(0, 1, n1)]
gs_label = gs/nS
centers_label = centers/um

fig = figure(1, figsize=(6,7))
gs=gs*1e9
ax1 = subplot(231)

for k in range(n1): semilogx(centers/um, i_thresholds[:,k], color=colors[k])

ax1.set_xlabel(r'$\Delta$ x1/2 ($\mu$m)') 
ax1.set_ylabel('Rheobase (pA)') 
ax1.set_ylim(50,72)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.xaxis.set_minor_formatter(ScalarFormatter())

ax2 = subplot(232)
for k in range(n1): semilogx(centers/um, v_thresholds[:,k], color=colors[k])

ax2.set_xlabel(r'$\Delta$ x1/2 ($\mu$m)') 
ax2.set_ylabel('$AP threshold$ (mV)') 
ax2.set_ylim(-59,-48)
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.xaxis.set_minor_formatter(ScalarFormatter())

ax3 = subplot(233)
for k in range(n1):  semilogx(centers/um, IR[:,k], color=colors[k])

ax3.set_xlabel(r'$\Delta$ x1/2 ($\mu$m)') #, fontsize=14)
ax3.set_ylabel('$Input resistance$ (MΩ)') 
ax3.set_ylim(270,350)
ax3.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax3.xaxis.set_minor_formatter(ScalarFormatter())

ax4 = subplot(234)
for k in range(n1): semilogx(centers/um, ff[:,k], color=colors[k])

ax4.set_xlabel(r'$\Delta$ x1/2 ($\mu$m)') #, fontsize=14)
ax4.set_ylabel('$Frequency$ (Hz)') 
ax4.set_ylim(-1,28)
ax4.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax4.xaxis.set_minor_formatter(ScalarFormatter())

ax5 = subplot(235)
for k in range(n1): semilogx(centers/um,np.abs(v_ahps[:,k]), color=colors[k], label=f'{test_param}={gs[k]:.2f} nS')
ax5.set_xlabel(r'$\Delta$ x1/2 ($\mu$m)') #, fontsize=14)
ax5.set_ylabel('$AHP amplitude $ (mV)') 
ax5.set_ylim(10,21)
ax5.legend(loc="upper left")
ax5.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax5.xaxis.set_minor_formatter(ScalarFormatter())

subplots_adjust(hspace=0.35, wspace=0.4)
# fig.savefig(f"AIS_with_distance_K_{test_param}_{gs}.pdf", bbox_inches='tight', dpi=300)
# fig.savefig(f"AIS_with_distance_K_{test_param}_{gs}.png", bbox_inches='tight', dpi=300)
show()    


