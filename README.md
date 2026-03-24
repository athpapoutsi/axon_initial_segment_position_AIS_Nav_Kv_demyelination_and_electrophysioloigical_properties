# Effect of axon initial segment (AIS) position, Nav and Kv AIs conductances and demyelination on electrophysioloigical properties
Code corresponding to the paper
"[Contactin2 Regulates Axonal Organization, Myelination, and Excitability of Hippocampal Somatostatin Interneurons]

Model is based on (https://elifesciences.org/articles/53432)" by Sarah Goethals and Romain Brette.
Adapted by Athanasia Papoutsi.

### Requirements

In addition to standard scientific packages (numpy, matplotlib), this code requires:
* joblib
* [Brian 2](http://briansimulator.org)


## Relation between excitability and AIS geometry / demyelination

* `AIS_location.py`: 
 We use a simple biophysical model of spike initiation in the AIS. With this model, we show how the electrophysioloigical properties depend on AIS position and demyelination.


## Robustness and modulation by AIS Nav and Kv condutances

* `AIS_location_paramters.py`: 
Using the same model, we show how AIS Nav and Kv conductances further modulate the electrophysiological features. 
