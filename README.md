# crosspy

Python framework for FFT-based image cross correlation. Uses subset method to obtain the deformation field.

Bergsmo & McAuliffe 2020

After Experimental Micromechanics group Matlab code ; eg Jiang et al - "Deformation compatability in a single crystalline Ni superalloy" - http://dx.doi.org/10.1098/rspa.2015.0690 

Fourier cross correlation after Guizar-Sicairos et al - "Efficient subpixel image registration algorithms" - https://doi.org/10.1364/OL.33.000156 

Uses pyFFTW framework - https://pyfftw.readthedocs.io/en/latest/ 

------

Go to crosspy/decks for example implementation.

We employ an object-oriented approach:
- Imset class contains metadata and loading methods for image data.
- DIC class can be instantiated from an Imset or an externally loaded numpy array. Cumulative or sequential displacement and strain calculation can then be performed.

Additionally:
- Subpixel Fourier registration functions are included in the .XCF module.
- The discrete Fourier transform employed can be user-specified. We recommend the 'pyfftw_numpy' argument.
- Rigid body translation, rotation and polynomial fitting correction algorithms
- Discontinuity tolerance is incorporated via a "Heaviside" implementation 

------

Install by command line navigating to /crosspy header directory, then run "pip install ." 
