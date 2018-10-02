Introduction
============

``Py4Incompact3D`` is a library for postprocessing data produced by Xcompact3D simulations.
The aim of this project is to facillitate automated postprocessing of Xcompact3D simulations by
providing, at first:

* Mesh class: this stores the domain data for the simulation
* Case class: this stores the information of the case: boundary conditions, fields etc.

With these building blocks, complex postprocessing tools may be built - for example, derivative
calculateors to compute the vorticity and Q-criterion given the velocity field.

It is hoped that users of Xcompact3D will find this library useful and contribute to its
development, for instance by adding additional functionality.
