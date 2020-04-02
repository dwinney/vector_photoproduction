# vector_photoproduction
Framework for vector meson photo-production amplitudes off a proton.

### Executables
Requires ROOT (tested with version 6.17).

To build any example script in the `executables` folder, for example `test.cpp`, use:

```
mkdir build
cd build
cmake ..
make test
````

##### psi_comparison
Comparison of the unpolarized cross sections for the photoproduction of the Psi(1s) and Psi(2s) states near threshold in the GlueX kinematics. Optional flags are `-c [0:180]` to change the center-of-mass frame angle (in degrees) and `-lab` whether to plot in terms of photon lab frame of center-of-mass frame energy.

### Amplitudes
The main object of interest is the abstract `amplitude` class. This allows you to build observables (so far only DXS) from helicity amplitudes. Implemented so far are those considered in [1,2]:

* Pomeron exchange amplitude (t-channel)
* Single baryon resonance (s-channel)

Incoherent (interfering) sums of amplitudes may be constructed through the `amplitude_sum` class (see for example `executables/tests/sum_test.cpp`).


### References
* [1] "Theoretical model of the phi meson photoproduction amplitudes" Lesniak and Szczepaniak [[arXiv:hep-ph/0304007]](https://arxiv.org/abs/hep-ph/0304007)
* [2] "Double Polarization Observables in Pentaquark Photoproduction" JPAC Collaboration [[arXiv:1907.09393]](https://arxiv.org/abs/1907.09393)
