#!/usr/bin/python3
import sys, os, re, getopt, numpy, random, argparse, itertools
import scipy, scipy.stats
import matplotlib.pyplot as plt
import subprocess 
import statistics 
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors

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

    algorithm = list(set([t[0] for t in data]))
    runs = sorted(list(set([int(t[1]) for t in data if int(t[1]) >= limit])))

    ax = plt.subplot(1,1,1)
    for alg in sorted(algorithm):
        means = [float(t[2])/int(t[1]) * 1e9 for t in data if t[0] == alg and int(t[1]) >= limit]
        stdev = [float(t[3])/int(t[1]) * 1e9 for t in data if t[0] == alg and int(t[1]) >= limit]
        ax.errorbar(runs, means, stdev, marker=next(marker), linestyle=next(linestyle), fillstyle=next(fillstyle), color=next(color), label=r'$' + alg + '$')

    handles, labels = ax.get_legend_handles_labels()
    print(labels)
    new_handles = [handles[0], handles[3], handles[2], handles[1]]
    new_labels = [labels[0], labels[3], labels[2], labels[1]]
    plt.legend(new_handles, new_labels, loc=0, ncol=1)
    plt.grid(True)
    plt.title(r"Comparison of sequential sampling algorithms")
    # plt.title(r"Running time for $N=2^{50}$ averaged over $2/n$ repetitions")
    plt.xscale("log", basex=2)
    plt.yscale("log")
    plt.ylim(ymax=10**2.3, ymin=0.0)
    plt.xlabel("Sample size $n$")
    plt.ylabel("Time per sample (ns)")
    plt.yticks([3, 5, 10, 20, 40, 80, 160])
    plt.gca().get_yaxis().set_major_formatter(plt.ScalarFormatter())
    plt.tight_layout()
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
print(limit)

# Init plots
marker = itertools.cycle(('^', 'x', 'd', 'o', 'o', '*', 's')) 
fillstyle = itertools.cycle(('full', 'full', 'none', 'full', 'none', 'full', 'full')) 
# linestyle = itertools.cycle(("-", "--", "-.", ":")) 
linestyle = itertools.cycle(("-")) 
color = itertools.cycle(map(lambda x : colors.cnames[x], ["red", "darkmagenta", "blue", "orange", "green", "darkblue", "darkgreen"])) 
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
plt.rc('font', serif='cm')
plt.rc('font', size=13)
pp = PdfPages(output_file + '.pdf')

# Title
print("##########################")
print(bcolors.HEADER + "Sampling (Direct Comparison)" + bcolors.ENDC)
print("##########################")

running_times = []

with open(input_file) as infile:
    for line in infile:
        if (line.strip().split()[0] != "#"):
            running_times.append(line.strip().split())

# Create plot
plot_scaling(running_times)

# Finalize
pp.close()
