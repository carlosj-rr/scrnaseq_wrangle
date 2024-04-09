import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
from itertools import combinations
import seaborn as sns

files_list=glob.glob("*nsENOGs.list")
keys=[ x.split("_")[0] for x in files_list]

d = {}
for i in range(len(keys)):
    key = keys[i]
    infile = files_list[i]
    d[key] = np.genfromtxt(infile,dtype=str)
# d now holds the species 2-letter code as the key, and its list of nsENOGs as its values

# Now we can organize the data for the bar chart. First, we get the total list of ENOGs and deduplicate it
l=[] # a list of lists with all the nsENOGs
for i in d:
        l.append(list(d[i]))

# flatten and deduplicate the list
ddl=np.unique([item for sublist in l for item in sublist])
# Now we define a function that returns, for an ENOG, a 1 if it is present in the corresponding species of 'd', or 0 otherwise:
def where_is_enog(enog_id, dix):
        l=[]
        for i in dix:
                l.append(int(enog_id in dix[i]))
        return(l)
# This function can now be applied across all ENOGs to see where each one is being expressed - a presence/absence table, based on the lists from the dictionary
x = np.array([ where_is_enog(i, d) for i in ddl ])
data_table=pd.DataFrame(x, index=ddl, columns=d.keys()) # make the table
data_table['row_sum'] = data_table.sum(axis=1) # add row sums as a new column
data_table.sort_values('row_sum', ascending=False, inplace=True) # sort so as to put the most represented genes on top
data_table.to_csv("Out_Binary_table.csv",sep=",") # save the table as a csv
reduced=data_table.iloc[0:100,0:8] # create a reduced version for plotting
# Make a heatmap using seaborn
colnames = list(reduced.columns)
rownames = list(reduced.index)
as_np = reduced.to_numpy()
sns.set(font_scale=0.3)
ax = sns.heatmap(as_np, linewidths=0.01, annot=as_np, xticklabels=colnames, yticklabels=rownames, cmap="crest")
ax.xaxis.tick_top()
plt.savefig("Heatmap.svg")
plt.close()
"""
The next one is a neat function that returns the ENOGs shared beteween an arbitrary amount of
species ids from the dictionary. It does this by subsampling the 'data_table' object to include only
the species requested, then outputting the indices that corresponds to rows in which the sum of
the row was equal to the number of species.
"""
d_groups = {'Eumetazoan': ['Dr', 'Dm', 'Ml', 'Mm', 'Nv', 'Hv', 'Sp', 'Cg'],
 'Cnidarians': ['Hv', 'Nv', 'Sp'],
 'Coelenterates': ['Hv', 'Ml', 'Nv', 'Sp'],
 'Bilaterians': ['Dr', 'Dm', 'Mm', 'Cg'],
 'Verts+cnids': ['Dr', 'Hv', 'Mm', 'Nv', 'Sp'],
 'Bilats+cteno': ['Dr', 'Dm', 'Ml', 'Mm', 'Cg'],
 'Protostomes': ['Cg', 'Dm'],
 'Vertebrates': ['Dr', 'Mm'],
 'Verts+cteno': ['Dr', 'Ml', 'Mm']}

def shared_enogs(species: list, data_table: pd.DataFrame):
        tot = len(species)
        shared_enogs=np.array(data_table.index[data_table[species].sum(axis=1) == tot])
        return(shared_enogs)

# Now I want to do this for all possible combinations of species
keys_set=set(d.keys())
combi_list=[]
for n in range(2,len(keys_set)+1): # Starts with two, because ofc one species will share all its ENOGs with itself.
        curr_comb=list(combinations(keys_set,n))
        for i in curr_comb:
                combi_list.append(list(i))
idcs=[ "_".join(i) for i in combi_list ] # making the row indices
vals=list(map((lambda x: shared_enogs(x, data_table).size),combi_list)) # making the values by applying 'shared_enogs()' over all possible species combinations
combis_enog_counts=pd.DataFrame(vals, index=idcs) # chuck results into a pandas DF
combis_enog_counts.columns=["shared ENOGs"]
reordered=combis_enog_counts.sort_values("shared ENOGs",ascending=False) # sort in ascending order for a neater graph

width=0.5
ax = reordered.plot(kind="bar", title='Shared nsENOGs per species combinations', xlabel='Organism set', xticks=[], ylabel='count shared ENOGs', width=width, figsize=(8,8))
#labels=[ x.replace("_",", ") for x in reordered.index] # version that uses the list of organism IDs
# this other option adds the observed counts to the bar labels.
labels = [ str(reordered.index[i])+" ("+str(reordered.iloc[i,0])+")" for i in range(len(reordered.index)) ]

def autolabel(labels,width):
        for i in range(len(labels)):
                ax.text(i-(width/2), reordered.iloc[i,0], labels[i], ha='left', va='bottom', rotation=80., fontsize='small')
autolabel(labels, width)

plt.savefig("Bar_graph.png")
plt.close()
