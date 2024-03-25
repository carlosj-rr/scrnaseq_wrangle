import numpy as np
import matplotlib.pyplot as plt
import argparse

"""
Takes in for two species:
1. a list of an arbitrary system's eggNOG orthogroups (expENOGs - nervous system in this project
except for testing dissimilar systems)
2. a list of their total eggNOG orthogroups (totENOGs)

Then it does 100 permutations to estimate more or less how many ENOGs would be found to be
shared in both genomes by pure chance, given a sample the size of expENOGs and plots that as a
histogram (scaled down by the smallest observed expENOGs set in the study).
It plots the observed value of actual observed overlap in between both expENOGs, and the totENOGs.
"""
parser = argparse.ArgumentParser()
parser.add_argument("-exp_enogs1","--expENOGs_1", help="List of eggNOG orthogroups found in the nervous system of a species")
parser.add_argument("-exp_enogs2","--expENOGs_2", help="List of eggNOG orthogroups found in the entire genome of the same species")
parser.add_argument("-tot_enogs1","--totENOGs_1", help="Proportion of the full datasets that will be subsampled, number between 0 and 1")
parser.add_argument("-tot_enogs2", "--totENOGs_2", help="Filename of output figure. Extension will determine format (.png, .jpg, .svg, etc. - uses matplotlib.pyplot)")
parser.add_argument("-n_perms", "--number_of_permutations", help="Number of permutations for the histogram data")
parser.add_argument("-o", "--output_image_file", help="Filename for the output image - uses matplotlib.pyplot, so .png, .jpg, and .svg are accepted", default="Out.png")
parser.add_argument("-sc", "--scale_count", help="Scale to be applied to calculations so that all histograms are aligned - usually the lowest number of expENOGs in the whole study", default=1)
args=parser.parse_args()

exp_enogs1 = args.expENOGs_1
exp_enogs2 = args.expENOGs_2
tot_enogs1= args.totENOGs_1
tot_enogs2 = args.totENOGs_2
n_perms = args.number_of_permutations
output_file = args.output_image_file
scale_count = int(args.scale_count)


def compare_shared(set1,set2,curr_scale):
        shared=np.round(np.intersect1d(set1,set2).size/set1.size*curr_scale*100,1)
        return(shared)

def one_permutation(ns_a, ns_b, tot_a,tot_b, curr_scale):
        rng = np.random.default_rng()
        perm_a = rng.choice(tot_a,size=ns_a.size, replace=False)
        perm_b = rng.choice(tot_b,size=ns_b.size, replace=False)
        return(compare_shared(perm_a,perm_b, curr_scale))

#scale_count = 1826 # Nematostella's count of nsENOGs (the one with the least nsENOGs)

prefix1 = exp_enogs1.split("_")[0]
ns1 = np.genfromtxt(exp_enogs1,dtype=str)
tot1 = np.genfromtxt(tot_enogs1,dtype=str)

prefix2 = exp_enogs2.split("_")[0]
ns2 = np.genfromtxt(exp_enogs2,dtype=str)
curr_scale = scale_count/ns2.size
tot2 = np.genfromtxt(tot_enogs2,dtype=str)
obs_exp = compare_shared(ns1, ns2, curr_scale)
obs_tot = compare_shared(tot1, tot2, curr_scale)
perms_vals = np.array([ one_permutation(ns1, ns2, tot1, tot2, curr_scale) for i in range(100) ])

ax = plt.gca()
xlim_val = np.mean(perms_vals)*2
ax.set_xlim([0,xlim_val])
plt.hist(perms_vals, bins=10)
#plt.axvline(obs_exp, color="r")
plt.axvline(obs_tot, color="k")
plt.title("Expected % of shared OGs between "+ prefix1 + " and " + prefix2 + " (blue, N=100)\nwith observed % of shared nervous system OGs (red)", loc='left', fontsize='small')
plt.xlabel("% shared OGs")
plt.ylabel("Counts")
plt.savefig(output_file)
plt.close()