# vector_photoproduction
Framework for building amplitudes involving vector meson production via quasi-elastic scattering of a real photon on a nucleon target.

<p align="center">
  <img width="300" src="./doc/FeynmanDiagram.png">
</p>

Such processes are of interest at many experiments at JLab and the future EIC.

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

Observables are evaluated in terms of the center-of-mass energy, _s ,_ and cosine of the scattering angle, _zs_. Alternatively to easily interface with event generators, Lorentz vectors may be passed using the `event` struct. For example:
```c++
// In s-channel variables
double dxs = amp.differential_xsection(s, zs);

// or four-vectors
TLorentzVector pGamma, pTarget, pVector, pRecoil;
event fvecs(pGamma, pTarget, pVec, pRecoil);
double dxs = amp.differential_xsection(fvecs);
```

## EXECUTABLES
Requires [ROOT](https://root.cern.ch/) (tested with version 6.17) with [*MathMore*](https://root.cern.ch/mathmore-library) libraries installed.

To build any example script in the `executables` folder, for example `test.cpp`, use:

```bash
mkdir build
cd build
cmake ..
make test
````

All executables have the following two optional flags for customizing the plotted output
```bash
-f string             # Desired filename of output (default: executable_name.pdf)
-y [double:double]    # Manually set the y-range for plotting
```

#### [polarized_pentaquark](./executables/polarized_pentaquark.cpp)
Sensitivity study of double polarized observables to the LHCb pentaquarks in Hall A at JLab.
Reproduces the results in [2]. See [JPAC page on γp→J/ψp](http://cgl.soic.indiana.edu/jpac/polarizedPenta.php).

#### [asymmetry_pentaquark](./executables/asymmetry_pentaquark.cpp)
Sensitivity study of the beam asymmetry to the LHCb pentaquarks at GlueX at JLab.

Optional command line flags are:
```bash
-e double     # Fixed CM energy in GeV (default: 4.45 GeV)
-lab          # Toggle whether input energy value (see above) is in lab frame (default: false)
-10q          # Toggle plotting two P_c states at the same BR (default: false)
```
Outputs beam asymmetry for single P_c(4450) state at fixed energy as a function of theta or if `-10q` is passed, plots individual signal / background components and incoherent sum as a function of CM angle for P_c(4450) and P_c(4380) at equal BR.

#### [psi_comparison](./executables/psi_comparison.cpp)
Comparison of the unpolarized cross sections for the photoproduction of the Psi(1S) and Psi(2S) states near threshold in the GlueX kinematics.

Optional command line flags are:
```bash
-c double    # CM scattering angle in degrees (default: 0)  
-lab          # Toggle plotting with lab photon energy on x-axis (default: false)
```
Outputs PSI(2S) unpolarized cross-section and ratio of (2S)/(1S) plotted in pdfs.

#### [chi_c1_photoproduction](./executables/chi_c1_photoproduction.cpp)
Analytical model for the unpolarized cross section near threshold of axial vector states. Decomposed into different exchanges in the t-channel (e.g. omega, rho, phi).

Optional command line flags are:
```bash
-c double    # CM scattering angle in degrees (default: 0)
-integ        # Toggle integrated xsection instead of differential
```
Outputs the differential (or integrated) cross-section to a pdf.

#### [X3872_photoproduction](./executables/X3872_photoproduction.cpp)
Prediction for the unpolarized cross-section for exclusive X(3872) photoproduction at low momentum transfer and high energies of interest for the future EIC.

Optional command line flags are:
```bash
-c double    # CM scattering angle in degrees (default: 0)
-n int      # Number of points to plot with (default: 100)
-integ        # Toggle integrated xsection instead of differential
```

#### [Zc_photoproduction](./executables/Zc_photoproduction.cpp)
Unpolarized cross-section predictions for the photoproduction of the Z_c+(4200) by charged pion exchange. Reproduces the results in [3].

Optional command line flags are:
```bash
-n int        # Number of points to plot with (default: 100)
-m double     # max energy W to plot (default: 25)
-diff         # Toggle differential xsection instead of integrated
```

## PLOTTING
Plots are automatically created using the JPAC collaboration style guidelines. For more information see the [jpacStyle](https://github.com/dwinney/jpacStyle) library.

<p align="center">
  <img width="275" src="./doc/JPAClogo.png">
</p>

## REFERENCES
* [1] "Theoretical model of the phi meson photoproduction amplitudes" Lesniak and Szczepaniak [[arXiv:hep-ph/0304007]](https://arxiv.org/abs/hep-ph/0304007)
* [2] "Double Polarization Observables in Pentaquark Photoproduction" JPAC Collaboration [[arXiv:1907.09393]](https://arxiv.org/abs/1907.09393)
* [3] "Photoproduction of the charged charmoniumlike Z+c(4200)" Wang, Chen, and Guskov [[arXiv:1503.02125]](https://arxiv.org/abs/1503.02125)
