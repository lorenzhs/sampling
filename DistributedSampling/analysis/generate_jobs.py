#!/usr/bin/python
import sys, math, os, re, getopt, random, argparse

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
    output = job_folder + "/job_" + str(n) + "_" + str(N) + "_" + str(P)
    run_output = run_folder + "/run_" + str(n) + "_" + str(N) + "_" + str(P)
    out = open(output, 'w')
    out.write("# @ job_name = sample_" + str(n) + "_" + str(N) + "_" + str(P) + "\n") 
    out.write("# @ comment = 'test sample weak-scaling'" + "\n")
    out.write("# @ error = " + log_folder + "/$(job_name).$(jobid).out" + "\n")
    out.write("# @ output = " + log_folder + "/$(job_name).$(jobid).out" + "\n")
    out.write("# @ environment = COPY_ALL" + "\n")
    out.write("# @ wall_clock_limit = 00:10:00" + "\n")
    out.write("# @ notification = error" + "\n")
    out.write("# @ notify_user = seba.lamm@gmail.com" + "\n")
    out.write("# @ job_type = bluegene" + "\n")
    out.write("# @ bg_size = " + str(P/8) + "\n")
    out.write("# @ queue" + "\n")
    out.write("\n") 
    out.write("runjob --np " + str(P) + " --ranks-per-node 8 : ~/code/build/run_method" + method + " -n " + str(n) + " -N " + str(N) + " -k " + str(k) + " -i " + str(iterations) + " -output " + run_output)
    out.close()

#######################
# Main
#######################

parser = argparse.ArgumentParser(description='Generate jobs for parallel sampling.')
parser.add_argument('-P', type=int, nargs='+')
parser.add_argument('-n', type=int, nargs='+')
parser.add_argument('-N', type=int, nargs='?')
parser.add_argument('-k', type=int, nargs='?')
parser.add_argument('-i', type=int, nargs='?')
parser.add_argument('-jf', type=str, nargs='?')
parser.add_argument('-rf', type=str, nargs='?')
parser.add_argument('-lf', type=str, nargs='?')
args = parser.parse_args()

i = args.i
P = args.P
n = args.n
N = args.N
k = args.k
job_folder = args.jf
run_folder = args.rf
log_folder = args.lf

if not os.path.exists(job_folder):
    os.makedirs(job_folder)

if not os.path.exists(run_folder):
    os.makedirs(run_folder)

if not os.path.exists(log_folder):
    os.makedirs(log_folder)

# Title
print("##########################")
print(bcolors.HEADER + "Sampling (Generate Jobs)" + bcolors.ENDC)
print("##########################")

for sample_size in n:
	for p in P:
	    total_n = int(math.log(p, 2)) + sample_size
	    iterations = int((2**i)/(2**total_n)*p)
	    sample(total_n, N, k, p, iterations, "P")


