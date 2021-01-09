#`crosspy`

This `python` package performs efficient subpixel digital image correlation. Installation instructions can be found at the package header page [here](../README.md).

Two principal classes are used in this implementation: 
- `Imset`: holds information and meta-information for a stack of images.
- `DIC`: handles your DIC run, settings, and stores results.

Source code for these classes is found [here](../crosspy/Classes.py).

Example notebooks detailing example workflows can be found [here](../decks).

## Key classes, methods, attributes:

## `Imset`
This instantiates a DIC image class, holding metadata. 
The image is not loaded until the method .imload() is called.

#####Inputs:
- folder_path (pathlib.Path object) - path to the folder of images.
- extension (str) - image file extension to search for.
- indices (list or int) - specify which images to take.

#####Methods:
- self.imload(numbers) - Load and return corresponding images identiifed with metadata.

#####Attributes:
- self.folder - file path to saving folder.
- self.extension - extension of images.
- self.paths - file paths of images.
- self.names - file names of images.
- self.n_ims - number of images in the set.

### Imset.imload():
Loads the images identified in the parent Imset which are enumerated by the input.

#####Inputs
- numbers (list of int): Which images in the imset to load.

#####Outputs:
- imarray_stack: numpy array of images of dimensions [H, W, n_images]

## `DIC`
A class that holds information relating to a DIC run on an imageset
At the moment it runs on the top left portion of the image (ie. if there are pixels to the right and down 
that can't be accommodated by square subsets of given roi_size, they will be ignored).

#####Inputs:
- images (Imset) - Imset containing metadata for loading image stacks.
- roi (dict) - ROI settings for DIC. 
                This should be dictionary with keys "size_pass", "overlap_percentage", "XCF_mesh".
- filter_settings (list) - High / low pass Fourier filter settings. 
                This should be [high_pass_filter, high_pass_filter_width, low_pass_filter, low_pass_filter_width] in px.
- savingfolder (str or pathlib.Path) - location to write results to if requested.

#####Methods:
- self.run_sequential(cores, ffttype, hs, cc_t) - run DIC on images with displacements/strains relative to previous image in stack. Saves results to attributes.
- self.run_cumulative(cores, ffttype, hs, cc_t) - run DIC on images with displacements/strains relative to first image in stack. Saves results to attributes.

These can only be run after .run_sequential() or .run_cumulative()
- self.calculate_strain(strain_method) - infer the strain using calculated displacements. Saves results to attributes.
- self.correct(method, printing, fn) - correct the images by removing the calculated displacements. Returns an updated stack of images.
- self.plot_displacements(colmap) - plot the calculated displacements if available in attributes.

These can only be run after .calculate_strains()
- self.plot_strains(colmap) - plot the calculated strains if available in attributes.
- self.plot_strain_meta(bins) - plots histograms of strains if available in attributes.
- self.save_data(output_folder) - save results to hdf5.

#####Attributes:

Available after .run_sequential() or .run_cumulative()
- self.dx_maps - x displacement maps
- self.dy_maps - y displacement maps
- self.ph_maps - cross-correlation peak height maps

and if heaviside being used
- self.rd_maps
- self.th_maps
- self.hs_maps
- self.js_maps

Available after .calculate_strains()
- self.strain_11 - 11 (xx) strain
- self.strain_22 - 22 (yy) strain
- self.strain_12 - 12 (xy) strain
- self.strain_eff - von Mises effective strain
- self.rotation - 12 (xy) rotation
- self.deformation_gradient - full deformation gradient tensor


###DIC.run_sequential() or DIC.run_cumulative()
sequential: Run DIC with reference images as previous in the stack. Generates attributes for this class.
cumulative: Run DIC with reference images as first in the stack. Generates attributes for this class.

#####Inputs
- cores (int) - how many CPU cores to run subsets in parallel on.
- ffttype (str) - one of 'fftw_numpy', 'fftw_scipy', or 'numpy' - which FFT implementation to use.
- hs (bool) - run Heaviside or not.
- cc_t (float) - cross-correlation peak height threshold.

#####Outputs
- self.rd_maps
- self.th_maps
- self.hs_maps
- self.js_maps
- self.ph_maps
- self.dx_maps
- self.dy_maps

###DIC.correct()

Remove displacements saved to DIC attributes from the input Imset by fitting and subtracting a polynomial. Returns corrected images.

#####Inputs:
- method (str) - one of 'polynomial' or 'rigid'  displacement correction, latter is affine.
- printing (bool) - whether to print to console.
- fn - optional functional form to solve for BG correction. If none, quadratic assumed. See polynom_im_correct documentation for more.

#####Outputs:
- images_corrected - numpy array of corrected images. 



