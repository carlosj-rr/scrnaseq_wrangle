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
parser.add_argument("-sp1","--species_1", help="List of eggNOG orthogroups found in the nervous system of a species")
parser.add_argument("-sp2","--species_2", help="List of eggNOG orthogroups found in the entire genome of the same species")
#parser.add_argument("-tot_enogs1","--totENOGs_1", help="Proportion of the full datasets that will be subsampled, number between 0 and 1")
#parser.add_argument("-tot_enogs2", "--totENOGs_2", help="Filename of output figure. Extension will determine format (.png, .jpg, .svg, etc. - uses matplotlib.pyplot)")
parser.add_argument("-n_perms", "--number_of_permutations", help="Number of permutations for the histogram data", default=100)
parser.add_argument("-o", "--output_image_file", help="Filename for the output image - uses matplotlib.pyplot, so .png, .jpg, and .svg are accepted", default="Out.png")
parser.add_argument("-sc", "--scale_count", help="Scale to be applied to calculations so that all histograms are aligned - usually the lowest number of **expENOGs** in the whole study", default=1)
#parser.add_argument("-sc_tot", "--scale_count_tot", help="Scale to be applied to calculations so that all histograms are aligned - usually the lowest number of **totENOGs** in the whole study", default=1)
args=parser.parse_args()

sp1_name = args.species_1
sp2_name = args.species_2
sp1_prefix = sp1_name.split( )[0][0] + sp1_name.split( )[1][0]
sp2_prefix = sp2_name.split( )[0][0] + sp2_name.split( )[1][0]

exp_enogs1 = sp1_prefix+"_nsENOGs.list"
exp_enogs2 = sp2_prefix+"_nsENOGs.list"
tot_enogs1= sp1_prefix+"_totENOGs.list"
tot_enogs2 = sp2_prefix+"_totENOGs.list"
n_perms = int(args.number_of_permutations)
output_file = args.output_image_file
scale_count_ns = int(args.scale_count)
#scale_count_tot = int(args.scale_count_tot)


def compare_shared(set1,set2,curr_scale):
        shared=np.round(np.intersect1d(set1,set2).size/set1.size*curr_scale*100,1)
        return(shared)

def one_permutation(ns_a, ns_b, tot_a,tot_b, curr_scale):
        rng = np.random.default_rng()
        perm_a = rng.choice(tot_a,size=ns_a.size, replace=False)
        perm_b = rng.choice(tot_b,size=ns_b.size, replace=False)
        return(compare_shared(perm_a,perm_b, curr_scale))

#scale_count = 1826 # Nematostella's count of nsENOGs (the one with the least nsENOGs)

#prefix1 = exp_enogs1.split("_")[0]
ns1 = np.genfromtxt(exp_enogs1,dtype=str)
tot1 = np.genfromtxt(tot_enogs1,dtype=str)

#prefix2 = exp_enogs2.split("_")[0]
ns2 = np.genfromtxt(exp_enogs2,dtype=str)
curr_scale_ns = scale_count_ns/ns2.size
tot2 = np.genfromtxt(tot_enogs2,dtype=str)
#curr_scale_tot = scale_count_tot/tot2.size
obs_exp = compare_shared(ns1, ns2, curr_scale_ns)
#obs_tot = compare_shared(tot1, tot2, curr_scale_tot)
perms_vals = np.array([ one_permutation(ns1, ns2, tot1, tot2, curr_scale_ns) for i in range(n_perms) ])
print(f"\n#### {sp1_name} VS {sp2_name} ####\nRandomized distribution Min = {np.min(perms_vals)}\nMean = {np.mean(perms_vals)}\nMax = {np.max(perms_vals)}\nObserved exp = {obs_exp}\n")
top_yval = np.histogram(perms_vals)[0].max()
xlim_val = max(np.mean(perms_vals)*2,100)
ax = plt.gca()
ax.set_xlim([0,xlim_val])
plt.hist(perms_vals, bins=10)
plt.axvline(obs_exp, color="r")
#plt.axvline(obs_tot, color="g")
plt.title("Expected % of shared OGs between "+ sp1_name + " and " + sp2_name + " (blue, N=" + str(n_perms) +")\nwith observed % of shared nervous system OGs (red)", loc='left', fontsize='small')
plt.text(obs_exp*1.01,top_yval*0.98,"Obs nsENOG ∩", fontsize='small')
plt.text(np.mean(perms_vals), top_yval*0.05, "Exp ∩",fontsize="small",color="orange",horizontalalignment="center")
plt.xlabel("% shared OGs")
plt.ylabel("Counts")
plt.savefig(output_file)
plt.close()

def guts_vs_ns_figmaker(ml_guts, ml_ns, ml_tot, sp2_ns, sp2_tot, sp2_name, sp2_prefix):
        obs_guts_ns = compare_shared(ml_guts, sp2_ns, 1)
        obs_ns_ns = compare_shared(ml_ns, sp2_ns, 1)
        permut_guts_ns = np.array([ one_permutation(ml_guts, sp2_ns, ml_tot, sp2_tot, 1) for i in range(100) ])
        permut_ns_ns = np.array([ one_permutation(ml_ns, sp2_ns, ml_tot, sp2_tot, 1) for i in range(100) ])
        mean_guts_ns = permut_guts_ns.mean()
        mean_ns_ns = permut_ns_ns.mean()
        top_yval = max(np.histogram(permut_guts_ns)[0].max(),np.histogram(permut_ns_ns)[0].max())
        xlim_val = min(mean_guts_ns*2, mean_ns_ns*2, 100)
        ax = plt.gca()
        ax.set_xlim([0,xlim_val])
        ax.set_ylim([0,top_yval*1.05])
        plt.hist(permut_guts_ns, bins=10, color="lightgreen")
        plt.axvline(obs_guts_ns, color="g")
        plt.hist(permut_ns_ns, bins=10, color="pink")
        plt.axvline(obs_ns_ns, color="r")
        plt.axvline((mean_guts_ns+mean_ns_ns)/2, color="k", ls="--")
        plt.title("Percent shared orthogroup sets between M. leidyi and "+sp2_name, loc="center", fontsize="medium")
        plt.xlabel("%shared OGs", fontsize="medium")
        plt.ylabel("Counts", fontsize="medium")
        plt.text(0, top_yval*0.5, "Mean green: " + str(mean_guts_ns) +"\nMean pink: " +str(mean_ns_ns), fontsize="small")
        plt.savefig("Mlguts_"+sp2_prefix+"ns.svg")
        plt.close()
        return