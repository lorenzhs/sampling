#!/usr/bin/python3
import sys, os, re, getopt, numpy, random, argparse
import scipy, scipy.stats
import matplotlib.pyplot as plt
import subprocess 
from matplotlib.backends.backend_pdf import PdfPages

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


#######################
# Main
#######################

parser = argparse.ArgumentParser(description='Evaluate different sampling methods.')
parser.add_argument('-i', type=str, nargs='?')
parser.add_argument('-f', type=str, nargs='?')
args = parser.parse_args()

input_file = args.i
output_file = args.f

# Title
print("##########################")
print(bcolors.HEADER + "Sampling (Direct Comparison)" + bcolors.ENDC)
print("##########################")

running_times = []

with open(input_file) as infile:
    for line in infile:
        args = re.split('=| ', line.strip())
        print(args)
        if len(args) > 10:
        # if len(args) > 10 and args[2] in ["mkl", "std"] and int(args[-1]) == 2**50:
            # running_times.append((args[2], int(args[-7]), float(args[4]) / 1000, float(args[6]) / 1000))
            running_times.append(("GPU", int(args[9]), float(args[1]) / 1000, float(args[3]) / 1000))
        # if (line.strip().split()[0] != "#"):
            # running_times.append(line.strip().split())

out = open(output_file, "w+")
out.write('\n'.join("%s %s %s %s" % x for x in sorted(running_times)))

