# Copyright 2018 Georgios Deskos, Paul Bartholomew
# Copyright 2021 University of Edinburgh
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed
# under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied. See the License for the
# specific language governing permissions and limitations under the License.

import Py4Incompact3D

from warnings import warn

import numpy as np

class Field():

    def __init__(self, *args, **kwargs):

        super().__init__()
        
        self.dtype = np.float64 # Default to double precision
        self.fromfile = True # By default load from a file

        if len(args) == 1:
            warn("You are using an old-style initialisation, the future is dynamic!", DeprecationWarning)

            instance_dictionary = args[0]
            
            self.name = instance_dictionary["name"]
            self.description = instance_dictionary["description"]

            properties = instance_dictionary["properties"]

            self.file_root = properties["filename"]
            self.direction = properties["direction"]

            if "precision" in properties:
                if properties["precision"] == "single":
                    self.dtype = np.float32
                else:
                    self.dtype = np.float64

            if "fromfile" in properties:
                self.fromfile = properties["fromfile"]
        else:
            for arg, val in kwargs.items():
                if arg == "name":
                    self.name = val
                elif arg == "description":
                    self.description = val
                elif arg == "file_root":
                    self.file_root = val
                elif arg == "direction":
                    self.direction = val
                elif arg == "precision":
                    if val == "single":
                        self.dtype = np.float32
                    else:
                        self.dtype = np.float64

        self.data = {}

    def _read(self, filename, nx, ny, nz, dtype=np.float64):
        """ Reads a datafile generated by Incompact3D into a (3D) numpy array. 
        
        :param filename: The file to read.
        :param nx: The mesh x resolution.
        :param ny: The mesh y resolution.
        :param nz: The mesh z resolution.
        :param dtype:
        """
        
        N = nx * ny * nz
        with open(filename, "rb") as bindat:
            fldat = np.fromfile(bindat, dtype)
            assert(len(fldat) == N)
            
        return np.reshape(fldat, (nx, ny, nz), "F")
    
    def _read_hdf5(self, t, nx, ny, nz):
        """ Reads a datafile generated by Incompact3D into a (3D) numpy array. 
        
        :param t:  The timestep to load.
        :param nx: The mesh x resolution.
        :param ny: The mesh y resolution.
        :param nz: The mesh z resolution.
        """
        import h5py # XXX: This doesn't seem right...
        
        N = nx * ny * nz
        with h5py.File(self.file_root, "r") as h5dat:
            h5path = "/Step" + str(t) + "/" + self.name
            arr = h5dat[h5path]
            fldat = np.zeros(arr.shape)
            arr.read_direct(fldat)
            fldat = fldat.flatten("F")
            assert(len(fldat) == N)

        return np.reshape(fldat, (nx, ny, nz), "F")
        
    def _to_fortran(self, time=-1):
        """ Converts data fields from internal (C) to Fortran ordering.
        
        :param time: The time(s) to convert.
        :type time: int or list of int
        """

        if time == -1:
            conv_times = self.data.keys()
        elif isinstance(time, int):
            conv_times = [time]
        elif isinstance(time, list):
            conv_times = time
        else:
            raise ValueError

        for t in conv_times:
            self.data[t] = np.swapaxes(self.data[t], 0, 2)

    def _from_fortran(self, time=-1):
        """ Converts data fields from Fortran to internal (C) ordering.

        :param time: The time(s) to convert.
        :type time: int or list of int
        """

        if time == -1:
            conv_times = self.data.keys()
        elif isinstance(time, int):
            conv_times = [time]
        elif isinstance(time, list):
            conv_times = time
        else:
            raise ValueError

        for t in conv_times:
            self.data[t] = np.swapaxes(self.data[t], 2, 0)

    def load(self, mesh, time=-1):
        """ Loads a datafield timeseries.

        :param mesh: The mesh the data is stored on.
        :param time: Time(s) to load data, default value -1 means load all times.

        :type mesh: Py4Incompact3D.postprocess.mesh.Mesh
        :type time: int or list of int
        """
        read_hdf5 = False

        # Find all files to load
        if time == -1:
            load_times = range(1000) # This corresponds to 4-digit timestamp
        elif isinstance(time, int):
            load_times = [time]
        elif isinstance(time, list):
            load_times = time
        else:
            raise ValueError

        if (Py4Incompact3D.HAVE_HDF5):
            import h5py # XXX: This doesn't seem right...
            if h5py.is_hdf5(self.file_root):
                read_hdf5 = True
            else:
                print(f"{self.file_root} is not an HDF5 file")

        if read_hdf5:
            for t in load_times:
                self.data[t] = self._read_hdf5(t, mesh.Nx, mesh.Ny, mesh.Nz)
        else:
            for t in load_times:
                zeros = ""
                read_success = False
                while (not read_success) and (len(zeros) < 10):
                    try:
                        filename = self.file_root + zeros + str(t)
                        self.data[t] = self._read(filename, mesh.Nx, mesh.Ny, mesh.Nz, self.dtype)
                    except FileNotFoundError:
                        msg = "Could not read " + filename + ", " + self.file_root + ", " + str(t)
                        print(msg)
                        zeros += "0"
                    else:
                        read_success = True

                if not read_success:
                    raise RuntimeError

    def _get_timestamp(self, t, timestamp_len=3):
        """ Set the timestamp for output according to format: phiXXX where XXX is the time
        left-padded with zeros.

        :param t: The time.
        :param timestep_len: The desired length of the timestamp.

        :type t: int
        :type timestep_len: int

        :returns: timestamp
        :rtype: str
        """

        timestamp = str(t)
        nzeros = timestamp_len - len(timestamp)
        if nzeros >= 0:
            timestamp = "0" * nzeros + timestamp
        else:
            msg = "Timestamp ({}) too short to format time {}".format(str(timestamp_len), str(t))
            raise RuntimeError(msg)

        return timestamp

    def write(self, time, timestamp_len=3):
        """ Output to binary file.

        :param time: The time(s) to write out to.
        :param timestep_len: How long should the timestep be? A length of 3 gives timestep 1 as 001,
        10 as 010 etc.

        :type time: int or list of int
        :type timestep_len: int
        """

        if time == -1:
            write_times = self.data.keys()
        elif isinstance(time, int):
            write_times = [time]
        elif isinstance(time, list):
            write_times = time
        else:
            raise ValueError

        for t in write_times:
            filename = self.file_root + self._get_timestamp(t, timestamp_len)

            # Dump to raw binary, numpy writes in 'C' order so we need to shuffle our
            # array so that 'C' order looks like 'Fortran' order...
            self._to_fortran(t)
            self.data[t].tofile(filename)

            # Shuffle back to 'C' order incase we want to keep working with the array
            self.data[t] = self._from_fortran(t)

    def clear(self):
        """ Cleanup data. """
        self.data= {}
        
