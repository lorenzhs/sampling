#!/usr/bin/python3
import sys, os, re, getopt, numpy, random, argparse, math
import scipy, scipy.stats
import matplotlib.pyplot as plt
import subprocess 
import statistics 
import itertools
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
# Plot data
#######################
def plot_scaling(data):
    print("##########################")
    print(bcolors.HEADER + "Generate Plot" + bcolors.ENDC)
    print("##########################")
    print("Generate plot...")

    bases = sorted(list(set([int(t[1]) for t in data if int(t[2]) == 1])))
    pes = sorted(list(set([int(t[2]) for t in data])))
    print(bases)
    print(pes)

    for base in bases:
        means = [float(t[3]) * 1e9 / (int(t[1])/int(t[2])) for t in data if int(t[1])/int(t[2]) == base]
        stdev = [float(t[4]) * 1e9 / (int(t[1])/int(t[2])) for t in data if int(t[1])/int(t[2]) == base]
        weak_scaling = means
        # weak_scaling = [float(means[0])/float(t) for t in means]
        # plt.errorbar(pes, weak_scaling, stdev, marker="^", label="$n/p$=" + str(base))
        plt.errorbar(pes, weak_scaling, marker=next(marker), linestyle='-', label=r'$n/P=2^{' + str(int(math.log(base, 2))) + '}$')

    plt.grid(True)
    plt.title(r"Weak scaling parallel sampling")
    # plt.title(r"Parallel efficiency for $N=2^{50}$ and $N/{np}$ repetitions")
    # plt.xlabel(r"$\log_{10}$(Sample size)")
    plt.xscale("log", basex=2)
    # plt.yscale("log")
    plt.xlabel("Number of PEs $P$")
    plt.ylabel(r"Running time / $(n/P)$ (ns)")
    plt.legend(loc=0)
    plt.tight_layout()
    # plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    pp.savefig()
    plt.close()


#######################
# Main
#######################

parser = argparse.ArgumentParser(description='Evaluate different sampling methods.')
parser.add_argument('-i', type=str, nargs='?')
parser.add_argument('-f', type=str, nargs='?')
parser.add_argument('-l', type=int, nargs='?')
args = parser.parse_args()

input_file = args.i
output_file = args.f
limit = args.l
# marker = itertools.cycle((".", ",", "o", "v", "^", "+", "x", "d")) 
marker = itertools.cycle(("o", "v", "x", "d")) 
linestyle = itertools.cycle(("-", "--", "-.", ":")) 
print(limit)

# Init plots
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
plt.rc('font', serif='cm')
plt.rc('font', size=13)
pp = PdfPages(output_file + '.pdf')

# Title
print("##########################")
print(bcolors.HEADER + "Sampling (Weak Scaling)" + bcolors.ENDC)
print("##########################")

running_times = []

with open(input_file) as infile:
    for line in infile:
        args = re.split('=| ', line.strip())
        print(args)
        running_times.append((args[2], args[4], args[8], args[10], args[12]))

# Create plot
plot_scaling(running_times)

# Finalize
pp.close()
