import numpy as np
import matplotlib.pyplot as plt
import argparse

"""
Takes in the set of nervous system eggNOG orthogroups (nsENOGs) and the same species' complete set
of ENOGs (totENOGs), plus a proportion to sample from them.

It excludes the nsENOGs from the totENOGs set, and randomly samples for each the proportion of the
dataset sample size established by the proportion passed as an argument. For each dataset
(ns and tot), it calculates the % overlap in orthogroups for two random subsets of the proportion.

It then outputs a plot with the distribution of the random values, and their means.
"""

parser = argparse.ArgumentParser()
parser.add_argument("-ns_enogs","--nsENOGs", help="List of eggNOG orthogroups found in the nervous system of a species")
parser.add_argument("-tot_enogs","--totENOGs", help="List of eggNOG orthogroups found in the entire genome of the same species")
parser.add_argument("-p","--sample_prop", help="Proportion of the full datasets that will be subsampled, number between 0 and 1")
parser.add_argument("-o", "--output_figname", help="Filename of output figure. Extension will determine format (.png, .jpg, .svg, etc. - uses matplotlib.pyplot)")

args=parser.parse_args()

ns_dataset = args.nsENOGs
tot_dataset = args.totENOGs
ssprop = float(args.sample_prop)
output_figfile = args.output_figname

def shared_percent(set1,set2):
        shared=np.round(np.intersect1d(set1,set2).size/set1.size*100,1)
        return(shared)

def self_permute(ns, tot, ssprop):
        ns_chunksize = np.round(ns.size*ssprop).astype(int)
        tot_chunksize = np.round(tot.size * ssprop).astype(int)
        rng = np.random.default_rng()
        ns_obs1 = rng.choice(ns, size=ns_chunksize, replace=False)
        ns_obs2 = rng.choice(ns, size=ns_chunksize, replace=False)
        tot_obs1 = rng.choice(tot, size=tot_chunksize, replace=False)
        tot_obs2 = rng.choice(tot, size=tot_chunksize, replace=False)
        prop_ns = shared_percent(ns_obs1, ns_obs2)
        prop_tot = shared_percent(tot_obs1, tot_obs2)
        return(prop_ns, prop_tot)
    
ns = np.genfromtxt(ns_dataset, dtype=str)
tot = np.genfromtxt(tot_dataset, dtype=str)
tot_non_ns = np.setdiff1d(tot, ns)

perms = np.array([ self_permute(ns, tot_non_ns, ssprop) for i in range(1000) ])
print(f"\nInput sample proportion = {ssprop}\n\nPermutation results:\n{perms}")
ns_mean, tot_mean = np.mean(perms, axis=0)
print(f"\n\nPermutation column-wise means:\n\nNC set = {ns_mean}\nEverything else = {tot_mean}")
plt.hist(perms[:,0], bins=20, color="green", alpha=0.75, lw=0.1) # green, 75% opaque
plt.hist(perms[:,1], bins=20, color="blue", alpha=0.75, lw=0.1) # blue, 75% opaque
plt.axvline(ns_mean, color="r") # red
plt.axvline(tot_mean, color="k") # black
ax = plt.gca()
topxlim = 10+(np.max(perms))
ax.set_xlim([0,topxlim])
plt.savefig(output_figfile)
print(f"Output image: {output_figfile}")
plt.close()