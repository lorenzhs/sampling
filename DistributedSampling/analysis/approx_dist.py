#!/usr/bin/python3
import sys, os, re, getopt, random, argparse, math
import subprocess 
import statistics 
import scipy, scipy.stats
import numpy as np
import matplotlib.pyplot as plt
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
# Degree Distribution
#######################
def plot_distribution():
    print("##########################")
    print(bcolors.HEADER + "Degree Distribution" + bcolors.ENDC)
    print("##########################")

    N_pow = 40
    n_pow = 20
    m_pow = N_pow-1 
    [N, n, m] = [2**N_pow, 2**n_pow, 2**m_pow]

    p = m/N
    q = 1 - p

    print("Generate hypergeometric distribution...")
    # print("Hypergeometric")
    # print("Mean: " + str(scipy.stats.hypergeom.mean(N, m, n)))
    # print("Median: " + str(scipy.stats.hypergeom.median(N, m, n)))
    # print("Var: " + str(scipy.stats.hypergeom.var(N, n, m)))
    rv = scipy.stats.hypergeom.rvs(N, m, n, size=1000000)
    plt.hist(rv, 1000, label="Hyp", alpha=0.3)

    # print("Binomial")
    # print("Mean: " + str(scipy.stats.binom.mean(n, p)))
    # print("Median: " + str(scipy.stats.binom.median(n, p)))
    # print("Var: " + str(scipy.stats.binom.var(n, p)))
    rv = scipy.stats.binom.rvs(n, p, size=1000000)
    plt.hist(rv, 1000, label="Bin", alpha=0.3)

    # print("Poisson")
    # print("Mean: " + str(scipy.stats.poisson.mean(n * p)))
    # print("Median: " + str(scipy.stats.poisson.median(n * p)))
    # print("Var: " + str(scipy.stats.poisson.var(n * p)))
    # rv = scipy.stats.poisson.rvs(n * (m/N), size=1000000)
    # plt.hist(rv, 1000, label="Pois", alpha=0.3)

    # print("Normal")
    # print("Mean: " + str(scipy.stats.poisson.mean(n * p)))
    # print("Median: " + str(scipy.stats.poisson.median(n * p)))
    # print("Var: " + str(scipy.stats.poisson.var(n * p)))
    # rv = scipy.stats.norm.rvs(n * p, scale=math.sqrt(n * p * q * (N-n)/(N-1)), size=1000000)
    # plt.hist(rv, 1000, label="Norm", alpha=0.3)

    plt.grid(True)
    plt.legend(loc=0)
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    pp.savefig()
    plt.close()


#######################
# Main
#######################

parser = argparse.ArgumentParser(description='Evaluate hypergeometric distribution.')
parser.add_argument('-f', type=str, nargs='?')
args = parser.parse_args()

out_file = args.f

# Init plots
pp = PdfPages(out_file + '.pdf')

# Network Generation
print("##########################")
print(bcolors.HEADER + "Hypergeometric Evaluation" + bcolors.ENDC)
print("##########################")

# Degree Distribution
plot_distribution()

# Finalize
pp.close()
