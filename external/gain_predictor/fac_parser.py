import numpy as np
import re
import warnings


class LevelsData:

    def get_level_index(self, level_key_string):
        for line in self.data:
            if str(line[-1]) == level_key_string:
                return line[0]

    def __init__(self, filename):
        self.data = np.genfromtxt(filename,
                                  skip_header=12,
                                  autostrip=True,
                                  dtype=None,
                                  )


class StructureData:
    def __init__(self, filename, levels_data, lower_level, upper_level):

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.data = np.genfromtxt(filename,
                                      skip_header=13,
                                      invalid_raise=False,
                                      autostrip=True,
                                      dtype=[
                                          ('upper_level_index', int),
                                          ('upper_level_j', int),
                                          ('lower_level_index', int),
                                          ('lower_level_j', int),
                                          ('transition_energy', float),
                                          ('gf', float),
                                          ('transition_rate', float),
                                          ('multipole', float)
                                      ],
                                      )

        lower_index = levels_data.get_level_index(lower_level)
        upper_index = levels_data.get_level_index(upper_level)
        for line in self.data:
            if line[0] == upper_index and line[2] == lower_index:
                self.data = line

    def __getitem__(self, item):
        return self.data[item]


class CollisionStrengthData:

    def find_grid_size(self):
        with open(self.filename) as file:
            match = re.search(r"NUSR\s=\s(\d*)", file.read())
            return int(match.group(1))

    def parse_line(self, lower_index, upper_index):
        with open(self.filename) as file:
            for i, line in enumerate(file):
                match = re.search(r"\s{{5}}{}\s{{2}}\d\s*{}\s{{2}}\d\s{{2}}(\d*.\d*E.{{3}})"
                                  .format(lower_index, upper_index), line)
                if match:
                    return float(match.group(1)), i
            raise Exception("Transition collisional strengths not found!")

    def parse_data(self):
        data = np.genfromtxt(self.filename,
                             skip_header=self.line_index + 2,
                             dtype=[('energy', float), ('collision_strength', float), ('cross_section', float)],
                             max_rows=self.grid_size)
        return data

    def __init__(self, filename, levels_data, lower_level, upper_level):
        self.filename = filename
        self.grid_size = self.find_grid_size()

        lower_index = levels_data.get_level_index(lower_level)
        upper_index = levels_data.get_level_index(upper_level)
        self.transition_energy, self.line_index = self.parse_line(lower_index, upper_index)
        self.data = self.parse_data()

    def __str__(self):
        return "Transition energy: " + \
               str(self.transition_energy) + "\n" + \
               "Energy[eV] CollStr CrossSec\n" + \
               str(self.data)

    def __getitem__(self, item):
        return self.data[item]
