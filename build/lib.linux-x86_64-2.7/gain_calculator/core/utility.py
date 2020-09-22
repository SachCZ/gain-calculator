# -*- coding: utf-8 -*-

import contextlib
import os
import sys


@contextlib.contextmanager
def no_stdout(output=os.devnull):
    sys.stdout.flush()  # <--- important when redirecting to files
    save_stdout = sys.stdout
    new_stdout = os.dup(1)
    devnull = os.open(output, os.O_WRONLY)
    os.dup2(devnull, 1)
    os.close(devnull)
    sys.stdout = os.fdopen(new_stdout, 'w')
    yield
    sys.stdout = save_stdout


# taken from https://gist.github.com/aubricus/f91fb55dc6ba5557fbab06119420dd6a
def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    """
    A helper function to be called in a loop to create terminal progress bar. Could be used for example to display
    progress of population calculation.

    :param int iteration: current iteration
    :param int total: total iterations
    :param str prefix: prefix string, '' by default
    :param str suffix: suffix string
    :param int decimals: positive number of decimals in percent complete
    :param int bar_length: character length of bar
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '█' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()
