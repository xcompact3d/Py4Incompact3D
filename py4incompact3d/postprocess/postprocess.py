# Copyright 2018 Georgios Deskos, Paul Bartholomew

# Licensed under the Apache License, Version 2.0 (the "License"); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software distributed
# under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied. See the License for the
# specific language governing permissions and limitations under the License.

from warnings import warn

from .input_reader import InputReader
from .mesh import Mesh
from .fields import Field

MESH_PROPERTIES=["n","l","bc","beta","stretched","yp"]

class Postprocess():

    """
    Postprocess is the highest level class of the Py4Incompact3D package. Import this class
    and instantiate it with a path to an input file to begin running Py4Incompact3D. Use
    the ``fields'' attribute to access other objects within the model.

    inputs:
        input_file: str - path to the nml input file
    outputs:
        self: post - an instantiated post object
    """

    def __init__(self, *args, **kwargs):

        if len(args) == 1:
            warn("You are using an old-style initialisation, the future is dynamic!", DeprecationWarning)

            input_file = args[0]
        
            self.input_reader=InputReader()
            self.input_file = input_file
            self.fields, self.mesh = self._process_input()

        else:
            self.mesh = Mesh(*args, **kwargs)
            self.fields = {}

    def add_field(self, name, filepath, **kwargs):

        description = ""
        direction = -1   # Assume field is a scalar
        dtype = "double" # Xcompact3d uses double by default
        io_name = "solution-io"
        
        self.fields[name] = Field(name=name, file_root=filepath,
                                  description=description,
                                  direction=direction,
                                  dtype=dtype,
                                  io_name=io_name)
        
    def _process_input(self):
        return self.input_reader.read(self.input_file)

    def load(self, **kwargs):
        """ Load data.
        """

        load_vars = self.fields.keys()
        time = -1
        for arg, val in kwargs.items():
            if "vars" == arg:
                load_vars = val
            elif "time" == arg:
                time = val

        for var in load_vars:
            if self.fields[var].fromfile:
                self.fields[var].load(self.mesh, time)

    def write(self, **kwargs):
        """ Write data.
        """

        vars = "all"
        time = -1
        for arg, val in kwargs.items():
            if "vars" == arg:
                vars = val
            elif "time" == arg:
                time = val

        if vars == "all":
            vars = self.fields.keys()
        print(vars)

        for var in vars:
            self.fields[var].write(time)

    def clear_data(self, vars="all"):
        """ Clear stored data fields. """

        if vars == "all":
            vars = self.fields.keys()

        for var in vars:
            self.fields[var].clear()

