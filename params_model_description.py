
"""
Parameters used for the simplified biophysical model.

"""

from brian2 import *

# Passive parameters
Cm = 0.9* uF / cm ** 2
gL = 1./(75000 * ohm * cm**2) 
Ri = 100. * ohm * cm

ENa = 70. * mV
EK = -90. * mV

# Channels kinetics
T = 33.
factor = (1. / 2.8) ** ((T - 23.) / 10.)

# Na channels parameters
Va_soma = -30.*mV 
Ka = 5.05 * mV
Taum_max = factor * 0.15 * ms

Vh_soma = -60.*mV
Kh = 5. * mV  
Tauh_max = factor * 2.* ms  

#For AIS
Va = -32.5 * mV  
Vh = -65.*mV   

# K channels parameters
Vn = -70.*mV 
Kn = 30.*mV 
Taun_max = 8. * ms 

# Na channels parameters
gna_soma = 1125. * (siemens / meter ** 2) 
gna_dend = 50. * (siemens / meter ** 2) 

# K channels parameters
gk_soma = 125. * (siemens / meter ** 2) 
gk_dend = 50. * (siemens / meter ** 2) 
