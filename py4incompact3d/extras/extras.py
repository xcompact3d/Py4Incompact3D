# Try to import h5py for use with HDF5 files
HDF5_AVAILABLE=False
try:
    import h5py
except ImportError:
    pass
else:
    HDF5_AVAILABLE = True

HAVE_DECOMP2D=False
try:
    import decomp2d
except ImportError:
    pass
else:
    HAVE_DECOMP2D = True

if HAVE_DECOMP2D:
    ADIOS2_AVAILABLE = bool(decomp2d.decomp4py.have_adios2)
else:
    ADIOS2_AVAILABLE = False

# Default to not using HDF5 or ADIOS2
HAVE_HDF5=False
HAVE_ADIOS2=False
    
