# Copyright 2018 Georgios Deskos

# Licensed under the Apache License, Version 2.0 (the "License"); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software distributed
# under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied. See the License for the
# specific language governing permissions and limitations under the License.

import numpy as np

class Field():

    def __init__(self, instance_dictionary):

        super().__init__()

        self.name = instance_dictionary["name"]
        self.description = instance_dictionary["description"]

        properties = instance_dictionary["properties"]
        self.file_root = properties["filename"]
        self.direction = properties["direction"]

    def load(self, time=None):
        """
        """

        if time == None:
            load_times = [""]
        elif isinstance(time, int):
            load_times = [time]
        elif isinstance(time, list):
            load_times = time
        else:
            raise ValueError

        data = {}
        for t in load_times:
            zeros = ""
            read_success = False
            while (not read_success) and (len(zeros) < 10):
                try:
                    filename = self.file_root + zeros + str(t)
                except:
                    zeros += "0"
                else:
                    read_success = True
                    
            if not read_success:
                raise RuntimeError

            with open(filename, "rb") as bindat:
                data[t] = np.fromfile(bindat, np.float64)

        return data
        
