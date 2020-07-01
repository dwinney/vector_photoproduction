# jpacPhoto
Framework for Object-Oriented Amplitude Analysis (OOAA) involving single meson production via quasi-elastic scattering of a real photon on a nucleon target. Focus on expandability and easy interfacing with Monte-Carlo tools and event generators.

<p align="center">
  <img width="300" src="./doc/FeynmanDiagram.png">
</p>

Such processes are of interest at many experiments at JLab and the future EIC.

Requires [ROOT](https://root.cern.ch/) (tested with version 6.17) with [*MathMore*](https://root.cern.ch/mathmore-library) libraries installed.

## AMPLITUDES
The main object of interest is the abstract [`amplitude`](./include/amplitudes/amplitude.hpp) class. This allows you to build [observables](./src/amplitudes/observables.cpp) from helicity amplitudes:

* Probability distribution ( √ Σ | A |^2 )
* Differential cross section ( dσ / dt )
* Integrated total cross section ( σ )
* Polarization asymmetries ( A_LL and K_LL )
* Spin density matrix elements ( ρ^α_λ,λ' )
* Integrated beam asymmetry ( Σ )
* Parity asymmetry ( P_σ )

Available amplitudes, so far, include:

* [Single baryon resonance](./include/amplitudes/baryon_resonance.hpp) (s-channel)
* [Pomeron exchange](./include/amplitudes/pomeron_exchange.hpp) (t-channel)
* [(fixed-spin and reggeized) Charged pseudo-scalar meson exchange](./include/amplitudes/vector_exchange.hpp) (t-channel)
* [(fixed-spin and reggeized) Vector meson exchange](./include/amplitudes/vector_exchange.hpp) (t-channel)
* [(fixed-spin) Dirac fermion exchange](./include/amplitudes/dirac_exchange.hpp) (u-channel)
* [(fixed-spin) Rarita-Schwinger fermion exchange](./include/amplitudes/rarita_exchange.hpp) (u-channel)

Incoherent (interfering) sums of amplitudes may be constructed through the [`amplitude_sum`](./include/amplitudes/amplitude_sum.hpp) class.

Observables are evaluated in terms of the invariant center-of-mass energy, s, and momentum transfer, t. Alternatively to easily interface with event generators, Lorentz vectors may be passed using the `event` struct to calculate s and/or t. For example:
```c++
// From four-vectors:
TLorentzVector pGamma, pTarget, pVector, pRecoil;
event fvecs(pGamma, pTarget, pVec, pRecoil);
double s = fvecs.s_man(), t = fvecs.t_man();

// Observable from s and t
double dxs = amplitude.differential_xsection(s, t);
```

## EXECUTABLES
To build any example script in the `executables` folder, for example `test.cpp`, use:

```bash
mkdir build
cd build
cmake ..
make test
````
Alternatively use `make JpacPhoto` to build a library file `libJpacPhoto.a` which may then be linked to other code to access header files and classes.

See documentation in each respective `.cpp` file for additional flag options.

#### [polarized_pentaquark](./executables/polarized_pentaquark.cpp)
Sensitivity study of double polarized observables, ALL and KLL, to the LHCb pentaquarks in Hall A at JLab.
Reproduces the results in [2]. See [JPAC page on γp→J/ψp](http://cgl.soic.indiana.edu/jpac/polarizedPenta.php).

#### [asymmetry_pentaquark](./executables/asymmetry_pentaquark.cpp)
Sensitivity study of the beam asymmetry to different LHCb pentaquark scenarios at GlueX at JLab.


#### [psi_comparison](./executables/psi_comparison.cpp)
Comparison of the unpolarized cross sections for the photoproduction of the Psi(1S) and Psi(2S) states near threshold in the GlueX kinematics.


#### [chi_c1_photoproduction](./executables/chi_c1_photoproduction.cpp)
Analytical model for the unpolarized cross section near threshold of axial vector states. Decomposed into different exchanges in the t-channel (e.g. omega, rho, phi).

#### [X3872_photoproduction](./executables/X3872_photoproduction.cpp)
Prediction for the unpolarized cross-section for exclusive X(3872) photoproduction at low momentum transfer and high energies of interest for the future EIC.

#### [Zc_photoproduction](./executables/Zc_photoproduction.cpp)
Unpolarized cross-section predictions for the photoproduction of the Z_c+(3900) and Z_c+(4200) by charged pion exchange. Reproduces the results in [3] with up to date widths and .

## PLOTTING
Plots are automatically created using the JPAC collaboration style guidelines. For more information see the [jpacStyle](https://github.com/dwinney/jpacStyle) library.

<p align="center">
  <img width="275" src="./doc/JPAClogo.png">
</p>

## REFERENCES
* [1] "Theoretical model of the phi meson photoproduction amplitudes" Lesniak and Szczepaniak [[arXiv:hep-ph/0304007]](https://arxiv.org/abs/hep-ph/0304007)
* [2] "Double Polarization Observables in Pentaquark Photoproduction" JPAC Collaboration [[arXiv:1907.09393]](https://arxiv.org/abs/1907.09393)
* [3] "Photoproduction of the charged charmoniumlike Z+c(4200)" Wang, Chen, and Guskov [[arXiv:1503.02125]](https://arxiv.org/abs/1503.02125)
