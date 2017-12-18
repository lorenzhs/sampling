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
    call = interface + " ../optimized/run_method" + method + " --exact_n --n=" + str(n) + " --N=" + str(N) + " --k=" + str(k) + " --seed=" + str(seed) + " --output=" + output
    print(call)

    subprocess.check_call(call, shell=True)
    with open(output) as infile:
        lines = infile.readlines()
        results[method].append(float(lines[0].split(" ")[2]))
    subprocess.check_call(["rm", output])


#######################
# Process results
#######################
def process_results(results, pes, data):
    for method in data:
        values = data[method]
        mean = statistics.mean(values)
        stdev = statistics.pstdev(values)
        results.append((method, pes, mean, stdev))

#######################
# Plot data
#######################
def plot_scaling(runs, data):
    print("##########################")
    print(bcolors.HEADER + "Generate Plot" + bcolors.ENDC)
    print("##########################")
    print("Generate plot...")

    p_means = [t[2] for t in data if t[0] == "P"]
    p_stdev = [t[3] for t in data if t[0] == "P"]
    weak_scaling = [float(p_means[0])/float(t) for t,e in zip(p_means, P)]

    plt.xticks(runs)
    plt.errorbar(runs, weak_scaling, marker="^", label="P (PE=" + str(P) + ")")

    plt.grid(True)
    plt.title(r"Parallel efficiency for $N=2^{" + str(N) + "}$ averaged over $" + str(seeds) + r"$ repetitions")
    # plt.xlabel(r"$\log_{10}$(Sample size)")
    plt.xlabel("Number of PEs")
    plt.ylabel("Weak scaling efficiency")
    plt.legend(loc=0)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    pp.savefig()
    plt.close()


#######################
# Main
#######################

parser = argparse.ArgumentParser(description='Evaluate weak scaling of parallel sampling.')
parser.add_argument('-f', type=str, nargs='?')
parser.add_argument('-P', type=int, nargs='+')
parser.add_argument('-n', type=int, nargs='?')
parser.add_argument('-N', type=int, nargs='?')
parser.add_argument('-k', type=int, nargs='?')
parser.add_argument('-s', type=int, nargs='?')
args = parser.parse_args()

output_file = args.f
seeds = args.s
P = args.P
n_per_pe = args.n
N = args.N
k = args.k

# Init plots
pp = PdfPages(output_file + '.pdf')

# Title
print("##########################")
print(bcolors.HEADER + "Sampling (Weak Scaling)" + bcolors.ENDC)
print("##########################")

# Output

# Perform test run for each seed
running_times = []
params = []

for p in P:
    times = { "P": [] }

    total_n = p * n_per_pe
    for seed in [random.randint(1, 1000) for num in range(seeds)]:
        sample(total_n, N, k, p, seed, "P", times)

    process_results(running_times, p, times)

# Create plot
plot_scaling(P, running_times)

# Finalize
pp.close()
