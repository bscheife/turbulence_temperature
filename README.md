# turbulence_temperature

This is one of the libraries I developed over the course of my PhD research to derive turbulent dissipation rates (ε and χ) from raw temperature microstructure measurements. There's a lot going on here, but the basic idea is to 
(i)   iteratively take chunks of the high resolution measurements (512 Hz), 
(ii)  transform them into fourier space, 
(iii) convert temperature spectra to temperature-gradient spectra
(iv)  integrate to obtain dissipation rates of temperature variance
(v)   find and subtract noise spectra from the measured signal
(vi)  fit theoretical spectra to the observed spectra using an MLE algorithm
(vii) estimate the Batchelor wavenumber from the fit, and from there derive the dissipation rate of turbulent kinetic energy

Check out https://doi.org/10.1029/2017JC013731 for details

## Getting started

The most important file, where most of the actual data processing happens, is in the function micT_eps_chi.m. This is called by its parent script t03_micT_processing.m. 

You can follow the basic workflow from start to finish by going through the files t01_xx.m, t02_xx.m, etc... in order.

## Acknowledgements

The basic structure of this code and workflow was provided by Barry Ruddick at Dalhousie University. The theory underlying the process is outlined in his paper *Maximum Likelihood Spectral Fitting: The Batchelor Spectrum*.

Jeff Carpenter at HZG in Germany iterated on this library before he passed it onto me. I fleshed out much of the code, made it more flexible or efficient, fixed bugs as necessary, and developed the documentation.
