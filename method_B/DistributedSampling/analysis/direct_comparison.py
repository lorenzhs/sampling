#!/usr/bin/python3
import sys, os, re, getopt, numpy, random, argparse
import scipy, scipy.stats
import matplotlib.pyplot as plt
import subprocess 
import statistics 
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
# Sampling procedures
#######################
def sample(n, N, k, P, seed, method, results):
    output = str(n) + "_" + str(N) + "_" + method + "_" + str(seed)
    if method in ["H", "D", "R"]:
        interface = "mpirun -n " + str(1) 
    elif method in ["P"]:
        interface = "mpirun -n " + str(P) 
    call = interface + " ../optimized/run_method" + method + " --n=" + str(n) + " --N=" + str(N) + " --k=" + str(k) + " --seed=" + str(seed) + " --output=" + output
    print(call)

    subprocess.check_call(call, shell=True)
    with open(output) as infile:
        lines = infile.readlines()
        results[method].append(float(lines[0].split(" ")[2]))
    subprocess.check_call(["rm", output])


#######################
# Process results
#######################
def process_results(results, sample_size, data):
    for method in data:
        values = data[method]
        mean = statistics.mean(values)
        stdev = statistics.pstdev(values)
        results.append((method, sample_size, mean, stdev))

#######################
# Plot data
#######################
def plot_scaling(runs, data):
    print("##########################")
    print(bcolors.HEADER + "Generate Plot" + bcolors.ENDC)
    print("##########################")
    print("Generate plot...")

    d_means = [t[2] for t in data if t[0] == "D"]
    d_stdev = [t[3] for t in data if t[0] == "D"]
    plt.errorbar(runs, d_means, d_stdev, marker="^", label="D")

    h_means = [t[2] for t in data if t[0] == "H"]
    h_stdev = [t[3] for t in data if t[0] == "H"]
    plt.errorbar(runs, h_means, h_stdev, marker="^", label="H")

    r_means = [t[2] for t in data if t[0] == "R"]
    r_stdev = [t[3] for t in data if t[0] == "R"]
    plt.errorbar(runs, r_means, r_stdev, marker="^", label="R")

    p_means = [t[2] for t in data if t[0] == "P"]
    p_stdev = [t[3] for t in data if t[0] == "P"]
    plt.errorbar(runs, p_means, p_stdev, marker="^", label="P (PE=" + str(P) + ")")

    plt.grid(True)
    plt.title(r"Running time for $N=2^{" + str(N) + "}$ averaged over $" + str(seeds) + r"$ repetitions")
    plt.xscale("log")
    # plt.xlabel(r"$\log_{10}$(Sample size)")
    plt.xlabel("Sample size")
    plt.ylabel("Time (s)")
    plt.legend(loc=0)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    pp.savefig()
    plt.close()


#######################
# Main
#######################

parser = argparse.ArgumentParser(description='Evaluate different sampling methods.')
parser.add_argument('-f', type=str, nargs='?')
parser.add_argument('-P', type=int, nargs='?')
parser.add_argument('-n', type=int, nargs='+')
parser.add_argument('-N', type=int, nargs='?')
parser.add_argument('-k', type=int, nargs='?')
parser.add_argument('-s', type=int, nargs='?')
args = parser.parse_args()

output_file = args.f
seeds = args.s
P = args.P
n = args.n
N = args.N
k = args.k

# Init plots
pp = PdfPages(output_file + '.pdf')

# Title
print("##########################")
print(bcolors.HEADER + "Sampling (Direct Comparison)" + bcolors.ENDC)
print("##########################")

# Output

# Perform test run for each seed
running_times = []
params = []

for sample_size in n:
    times = { "D": [], "H": [], "R": [], "P": [] }

    for seed in [random.randint(1, 1000) for num in range(seeds)]:
        sample(sample_size, N, k, P, seed, "D", times)
        sample(sample_size, N, k, P, seed, "H", times)
        sample(sample_size, N, k, P, seed, "R", times)
        sample(sample_size, N, k, P, seed, "P", times)

    process_results(running_times, sample_size, times)

# Create plot
plot_scaling([10**samples for samples in n], running_times)

# Finalize
pp.close()
