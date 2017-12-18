#!/usr/bin/python
import sys, os, re, getopt, random, argparse
import subprocess 
import statistics 

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
def sample(n, N, k, P, iterations, method, results):
    output = str(n) + "_" + str(N) + "_" + method + "_" + str(iterations)
    call = "../build/run_method" + method + " -n " + str(n) + " -N " + str(N) + " -k " + str(k) + " -i " + str(iterations) + " -output " + output
    print(call)

    subprocess.check_call(call, shell=True)
    with open(output) as infile:
        lines = infile.readlines()
        args = re.split('=| ', lines[0].strip())
        running_times.append((args[2], 2**n, float(args[4]), float(args[6])))
    subprocess.check_call(["rm", output])


#######################
# Write results
#######################
def write_results(runs, data):
    print("##########################")
    print(bcolors.HEADER + "Write results" + bcolors.ENDC)
    print("##########################")
    print("Generate plot...")

    out = open(output_file, "w+")

    d_means = [t[2] for t in data if t[0] == "D"]
    d_stdev = [t[3] for t in data if t[0] == "D"]
    for x in zip(runs, d_means, d_stdev):
        out.write("D %s %s %s\n" %x)

    h_means = [t[2] for t in data if t[0] == "H"]
    h_stdev = [t[3] for t in data if t[0] == "H"]
    for x in zip(runs, h_means, h_stdev):
        out.write("H %s %s %s\n" %x)

    r_means = [t[2] for t in data if t[0] == "R"]
    r_stdev = [t[3] for t in data if t[0] == "R"]
    for x in zip(runs, r_means, r_stdev):
        out.write("R %s %s %s\n" %x)

#######################
# Main
#######################

parser = argparse.ArgumentParser(description='Evaluate different sampling methods.')
parser.add_argument('-f', type=str, nargs='?')
parser.add_argument('-P', type=int, nargs='?')
parser.add_argument('-n', type=int, nargs='+')
parser.add_argument('-N', type=int, nargs='?')
parser.add_argument('-k', type=int, nargs='?')
parser.add_argument('-i', type=int, nargs='?')
args = parser.parse_args()

output_file = args.f
i = args.i
P = args.P
n = args.n
N = args.N
k = args.k

# Title
print("##########################")
print(bcolors.HEADER + "Sampling (Sequential Comparison)" + bcolors.ENDC)
print("##########################")

running_times = []
for sample_size in n:
    iterations = int((2**i)/(2**sample_size))
    sample(sample_size, N, k, P, iterations, "D", running_times)
    sample(sample_size, N, k, P, iterations, "H", running_times)
    sample(sample_size, N, k, P, iterations, "R", running_times)
    sample(sample_size, N, k, P, iterations, "SR", running_times)

write_results([2**samples for samples in n], running_times)

