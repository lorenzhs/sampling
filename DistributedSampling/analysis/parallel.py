#!/usr/bin/python
import sys, math, os, re, getopt, random, argparse
import subprocess 

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
def sample(n, N, k, P, iterations, method):
    output = str(n) + "_" + str(N) + "_" + method + "_" + str(P)
    interface = "job_submit -x+ -t 10 -T 10 -m 32000 -c p -p " + str(P) + "/1"
    call = interface + " mpirun ../build/run_method" + method + " -n " + str(n) + " -N " + str(N) + " -k " + str(k) + " -i " + str(iterations) + " -output " + output
    print(call)

    subprocess.check_call(call, shell=True)

#######################
# Main
#######################

parser = argparse.ArgumentParser(description='Evaluate weak scaling of parallel sampling.')
parser.add_argument('-f', type=str, nargs='?')
parser.add_argument('-P', type=int, nargs='+')
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
print(bcolors.HEADER + "Sampling (Weak Scaling)" + bcolors.ENDC)
print("##########################")

# Output

# Perform test run for each seed
running_times = []
params = []

for sample_size in n:
	for p in P:
	    total_n = int(math.log(p, 2)) + sample_size
	    iterations = int((2**i)/(2**total_n))
	    sample(total_n, N, k, p, iterations, "P")


