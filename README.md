# vector_photoproduction
Framework for vector meson photo-production amplitudes off a proton.

## Executables
Requires ROOT (tested with version 6.17).

To build any example script in the `executables` folder, for example `test.cpp`, use:

```
mkdir build
cd build
cmake ..
make test
````

##### polarized_pentaquark
Sensitivity study of double polarized observables to the LHCb pentaquarks in Hall A at JLab.
Reproduces the results in [2]. See [JPAC page on γp→J/ψp](http://cgl.soic.indiana.edu/jpac/polarizedPenta.php).

##### psi_comparison
Comparison of the unpolarized cross sections for the photoproduction of the Psi(1S) and Psi(2S) states near threshold in the GlueX kinematics. Optional flag are `-c [0:180]` to change the center-of-mass frame angle (in degrees) (default is forward scattering theta = 0).

##### chi_c1_photoproduction
Analytical model for the unpolarized cross section near threshold of axial vector states. Decomposed into different exchanges in the t-channel (e.g. omega, rho, phi).

## Amplitudes
The main object of interest is the abstract `amplitude` class. This allows you to build observables from helicity amplitudes:

* Differential cross section ( dσ / dt )
* Integrated total cross section ( σ )
* Polarization asymmetries ( A_LL and K_LL )

Available amplitudes, so far, are those considered in [1,2]:

* Single baryon resonance (s-channel)
* Pomeron exchange amplitude (t-channel)
* Vector meson exchange (t-channel)

Incoherent (interfering) sums of amplitudes may be constructed through the `amplitude_sum` class (see for example `executables/tests/sum_test.cpp`).


## References
* [1] "Theoretical model of the phi meson photoproduction amplitudes" Lesniak and Szczepaniak [[arXiv:hep-ph/0304007]](https://arxiv.org/abs/hep-ph/0304007)
* [2] "Double Polarization Observables in Pentaquark Photoproduction" JPAC Collaboration [[arXiv:1907.09393]](https://arxiv.org/abs/1907.09393)
