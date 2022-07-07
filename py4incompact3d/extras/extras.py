# Try to import h5py for use with HDF5 files
HAVE_HDF5=False
try:
    import h5py
except ImportError:
    pass
else:
    HAVE_HDF5 = True

HAVE_DECOMP2D=False
try:
    import decomp2d
except ImportError:
    pass
else:
    HAVE_DECOMP2D = True

if HAVE_DECOMP2D:
    HAVE_ADIOS2 = bool(decomp2d.decomp4py.have_adios2)
else:
    HAVE_ADIOS2 = False
