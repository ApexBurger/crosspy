# This is a document comprising desired updates to the code

Utility

    - Save and load functions for datasets -> Store in hdf5?
    - Folder GUI
    - Universal image type importer
    - Grayscale conversion for all image types
    - Filter GUI - port from matlab

Core functionality

    - Newton loop for re-mapping -> check (http://www.ncorr.com/index.php/dic-algorithms)
    - SEM correction algorithms

Performance

    - Parallelise strain calculation
    - On-the-go compilation of strain functions using numba - this requires a self written poly function -> numba speeds up these processes dramatically

Plotting functionality

    - EBSD overlays
    - Statistics