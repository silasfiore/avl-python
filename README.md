# python wrapper for Athena Vortex Lattice

This python wrapper for [AVL](https://web.mit.edu/drela/Public/web/avl/) can be used to read, write, and edit the text-based configuration files that define aircraft geometries. This can be especially useful for low-fidelity aerodynamic models of fixed-wing aircraft that need to be parametric. It could be used to automate the calculation of trim solutions or other flight performance data.

## avl.py

This script contains the main code with classes for each of the keywords that can be parsed by the fortran routines. Documentation on these input files can be found [here](https://web.mit.edu/drela/Public/web/avl/avl_doc.txt). The _splitspan_ method of the _surface_ class allows a lifting surface to be split into multiple segments. This can useful to insert control surfaces such as flaps, ailerons or elevons that span only a fraction of the total span.

## vsp2avl.py

This script can parse the .vsp xml format and convert it to an avl configuration. The VSP files contain aircraft geometries generated using the open source NASA program [OpenVSP](https://openvsp.org/). So far only lifting surfaces have been implemented. Openvsp can also be used to estimate the zero-lift drag of the aircraft. This drag coefficient can also be used as input for AVL to model the [aerodynamic drag](https://wright.grc.nasa.gov/airplane/drageq.html) of the aircraft.
