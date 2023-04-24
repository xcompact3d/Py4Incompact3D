Introduction
============

`Py4Incompact3D` is a library for postprocessig data produced by Xcompact3D simulations.
The aim of this project is to facillitate automated postprocessing of Xcompact3D simulations by
providing, at first:

* Mesh class: this stores the domain data for the simulation
* Case class: this stores the information of the case: boundary conditions, fields etc.

With these building blocks, complex postprocessing tools may be built - for example, derivative
calculateors to compute the vorticity and Q-criterion given the velocity field.

Installation
------------

* Clone the git repository
* Move into the location you cloned the code to and :code:`pip install .`, you may need to use
  :code:`pip install --user .` if you do not have write permissions in the default
  :code:`pip install` location.
* Test module can be imported by python interpreter: :code:`import py4incompact3d`
* [Optional] Install :code:`h5py` using :code:`pip3 install h5py` to work with HDF5 files generated
  by Xcompact3d - note this is only intended as a temporary solution until a proper integration with
  ADIOS2 is implemented.
  Upon import :code:`py4incompact3d` will report the status of HDF5 functionality, this can be toggled
  on/off by running `py4incompact3d.toggle_hdf5()` which will repeat the status report to confirm if
  HDF5 has been enabled.
  
Documentation
-------------

Documentation of functions can be found under `doc/build/latex/`.

To regenerate documentation, from the project root type :code:`make -C doc/ latexpdf` (requires
sphinx with the :code:`sphinx-fortran` extension `installation instructions`_).

.. _installation instructions: https://sphinx-fortran.readthedocs.io/en/latest/index.html

Contributing
------------

It is hoped that users of Xcompact3D will find this library useful and contribute to its
development, for instance by adding additional functionality.
