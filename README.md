# General
The Stochastic Pattern (i.e. random field) Generator (SPG) produces spatio-temporal pseudo-random Gaussian fields that satisfy the *proportionality of scales* property (Tsyroulnikov, Quart. J. Roy. Meteorol. Soc., 2001). The SPG is based on a third-order in time stochastic differential equation with a pseudo-differential spatial operator defined on a limited area 2D or 3D domain in the Cartesian coordinate system. The generated pseudo-random fields are homogeneous and isotropic in spacetime (with the scaled vertical and temporal coordinates). The correlation functions in any spatio-temporal direction belong to the Matern class. The spatio-temporal correlations are non-separable. A spectral-space numerical solver is implemented and accelerated exploiting properties of real-world geophysical fields, in particular, smoothness of their spatial spectra. The SPG is designed to create additive or multiplicative, or other spatio-temporal perturbations that represent uncertainties in numerical prediction models in geophysics.

References:

1) Tsyrulnikov M. and Gayfulin D. A Limited-Area Spatio-Temporal Stochastic Pattern Generator for ensemble prediction and ensemble data assimilation. Meteorol. Zeitschrift, 2016 (revised version under review).

2) Tsyrulnikov M. and Gayfulin D. A spatio-temporal stochastic pattern generator for simulation of uncertainties in geophysical ensemble prediction and ensemble data assimilation. ArXiv preprint, 2016, https://arxiv.org/abs/1605.02018.

3) Tsyrulnikov M. and Gayfulin D. A Stochastic Pattern Generator for ensemble applications. COSMO Technical Report N 29, (20 July) 2016, 51 pp, available at http://www.cosmo-model.org/content/model/documentation/techReports/docs/techReport29.pdf.

# Installation

Install using the standard cmake build procedure:

```bash
mkdir build && cd build
cmake ..
make && make check && make install
```

Or use [`nix`](https://nixos.org/nix/)!

# Usage

```bash
SPG.exe ./Config/testrun.cfg
```
