#!/usr/bin/env python

import matplotlib.pyplot as plt
import collections
import math
import sys


def parse_trace_file(path):
    """ 
    Simple spike file parsing function. Assumes no errors and might fail 
    silently on errors
    """
    trace = {"time":[], "voltage":[]}


    with open(path , "r") as f:


        for line in f:
            tokens = line.split(",")
            try:
                time = float(tokens[0].strip())
                voltage = float(tokens[1].strip())
            except:
                continue

            trace["time"].append(time)
            trace["voltage"].append(voltage)

    return trace

def main(path_trace=None):

    trace = parse_trace_file(path_trace)
    plot_trace(trace["time"], trace["voltage"])


def plot_trace(times, voltage):
    """
    Plot the histogram and target curve in one figure
    """
    plt.plot(times, voltage)

    plt.show()


if __name__ == "__main__":

    # Non checked command line parsing
    trace_to_parse = "./trace_0.0_vsoma.csv"
    if len(sys.argv) == 1:
        main(trace_to_parse)
    if len(sys.argv) == 2:
        main(sys.argv[1])
