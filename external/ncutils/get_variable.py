from netCDF4 import Dataset


def get_variable(filename, variable_name):
    root_group = Dataset(filename, "r", format="NETCDF4")
    data_dictionary = root_group.variables

    data = data_dictionary[variable_name][:][-1]

    root_group.close()

    return data
