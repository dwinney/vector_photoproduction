#   jpacPhoto
Framework for amplitude analysis involving single meson production via quasi-elastic scattering of a real photon on a nucleon target. Focus on expandability and easy interfacing with Monte-Carlo tools and event generators.

<p align="center">
  <img width="300" src="./doc/FeynmanDiagram.png">
</p>

Such processes are of interest at many experiments at JLab and the future EIC.

Requires [ROOT](https://root.cern.ch/) (tested with version 6.17) with [*MathMore*](https://root.cern.ch/mathmore-library) libraries installed.

##  INSTALLATION
To install clone normally and use:
```bash
mkdir build && cd build
cmake ..
cmake --build . --target install
```
This will create a `jpacPhoto/lib` with the linkable library. 


If you wish to also build the full suite of executables (e.g. to reproduce plots in [[1]](https://arxiv.org/abs/1907.09393) and [[2]](https://arxiv.org/abs/2008.01001)) you need to install the [jpacStyle](https://github.com/dwinney/jpacStyle) library and set environment variable as such:
```bash
# for bash
export JPACSTYLE=/path/to/jpacStyle

# for csh
setenv JPACSTYLE /path/to/jpacStyle
```

##  AMPLITUDES
The main object of interest is the abstract [`amplitude`](./include/amplitudes/amplitude.hpp) class. This allows you to build [observables](./src/amplitudes/observables.cpp) from helicity amplitudes:

* Probability distribution ( Σ_λ | A |^2 )
* Differential cross section ( dσ / dt )
* Integrated total cross section ( σ )
* Polarization asymmetries ( A_LL and K_LL )
* Spin density matrix elements ( ρ^α_λ,λ' )
* Integrated beam asymmetry ( Σ_4pi )
* Beam asymmetry in the y-direction ( Σ_y )
* Parity asymmetry ( P_σ )

Available amplitudes, so far, include:

* [Baryon resonance](./include/amplitudes/baryon_resonance.hpp) (s-channel)
* [Pomeron exchange](./include/amplitudes/pomeron_exchange.hpp) (t-channel)
* [(fixed-spin and reggeized) Charged pseudo-scalar meson exchange](./include/amplitudes/vector_exchange.hpp) (t-channel)
* [(fixed-spin and reggeized) Vector meson exchange](./include/amplitudes/vector_exchange.hpp) (t-channel)
* [Primakoff effect off nuclear target](./include/amplitudes/primakoff_effect.hpp) (t-channel)
* [(fixed-spin) Dirac fermion exchange](./include/amplitudes/dirac_exchange.hpp) (u-channel)
* [(fixed-spin) Rarita-Schwinger fermion exchange](./include/amplitudes/rarita_exchange.hpp) (u-channel)

Incoherent (interfering) sums of amplitudes may be constructed through the [`amplitude_sum`](./include/amplitudes/amplitude_sum.hpp) class.

##  REFERENCES
+ [1] [Double Polarization Observables in Pentaquark Photoproduction](https://arxiv.org/abs/1907.09393)
+ [2] [XYZ spectroscopy at electron-hadron facilities: Exclusive processes](https://arxiv.org/abs/2008.01001)
+ [3] [JPAC Website](http://cgl.soic.indiana.edu/jpac/index.php)

<p align="center">
  <img width="275" src="./doc/JPAClogo.png">
</p>