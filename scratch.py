import numpy as np
import pandas as pd

#read in the emapper annotations file
ann_dataset=pd.read_csv("Dr_refGenome-EGNG.emapper.annotations",sep="\t")
#read in the list of inferred nervous system genes
nc_genes=np.genfromtxt("../../Dr_ncGenes.csv",delimiter=",",dtype=str)

#see which one of the inferred NS genes have an annotation entry in the loaded dataset
ng_annotated=np.intersect1d(nc_genes,ann_dataset.iloc[:,0])

# make a list of all the row numbers in which a NS gene has its annotations
l=[]
for i in ng_annotated:
        exp=np.where(i == ann_dataset.iloc[:,0])[0]
        l.append(exp)

#Tidy list and make pretty:
l.sort()
l=np.array(l).flatten()

#BEAUTY!

np.any(np.in1d(nc_markers,ds_data.index))
np.any(np.in1d(nc_markers,ds_data.columns))

row_wise_sum=np.sum((dm_data > 0).astype(int),axis=1)
np.where(row_wise_sum >= dm_data.shape[1])[0].size


import nov_parser as novparser
np=novparser.np
pd=novparser.pd

ds_data=pd.read_csv("Mm_scRNAseq.csv")
nc_markers=np.array(pd.read_csv("NeuralMarkers_mouse.csv", sep=",")).flatten()

ds_data=novparser.correct_axes_orientation(ds_data,nc_markers)
nc_gene_ss=ds_data.loc[nc_markers,:]
gt0=np.sum((nc_gene_ss > 0).astype(int),axis=0)
novparser.top_tenperc_nc_gene_ids_extractor(gt0, ds_data, "Mm_ncGenes.csv")


import scipy as sc
import pandas as pd

ds_data=sc.io.mmread("matrix.mtx")
ds_data=ds_data.todense()
colnames=np.genfromtxt("barcodes.tsv",delimiter="\t",dtype=str)
for_rows=pd.read_csv("features.tsv",sep="\t",header=None)
rownames=np.array(for_rows.iloc[:,0])
cg_data=pd.DataFrame(ds_data,index=rownames,columns=colnames)

np.savetxt("Cg_bioMart_request.list",rownames,fmt='%s', delimiter=",")

cg_annots=pd.read_csv("eggNOGresults/Cg_refSeqs/Cg_refSeqs-EGNG.emapper.annotations",sep="\t")

cg_ann_summary=pd.DataFrame(list(zip(cg_annots.iloc[:,0],cg_annots.iloc[:,4])),columns=("gene_id","enogs"))

inferred_ncgenes=np.genfromtxt("Cg_ncMarkers_inferred.list",dtype=str)

cg_ann_summary.iloc[np.where(cg_ann_summary.loc[:,"gene_id"] == 'G33255')[0][0],1]

def enog_fetcher(enog_list,taxlev):
        for i in enog_list:
                if "@"+taxlev+"|" in i:
                        print("True")
                        hit=i
                        return(hit.split("@"+taxlev+"|")[0])
                else:
                        print("False")
                        return(None)
                
for i in enog_list:
        if "@"+taxlev+"|" in i:
             print("True")
             hit=i
        else:
             print("False")
             
np.setdiff1d(nc_markers,ds_data.index)
nc_markers=np.intersect1d(nc_markers,ds_data.index)

def enog_fetcher(enog_list,taxlev):
            for i in enog_list:
                if taxlev in i:
                        hit=i
                else:
                        hit=None
                if hit:
                        return(hit.split("@"+taxlev+"|")[0])
                else:
                        return(None)
                
                
emap_table=pd.read_csv("eggNOGresults/Cg_refSeqs/Cg_refSeqs-EGNG.emapper.annotations",sep="\t")
inferred_ncgenes=np.genfromtxt("Cg_ncMarkers_inferred.list",dtype=str)
in_common=np.intersect1d(inferred_ncgenes,cg_annots.loc[:,"#query"])
emap_table.index=emap_table.loc[:,"#query"]
emap_table.loc[in_common,"eggNOG_OGs"]
query_genes_enogs=emap_table.loc[in_common,"eggNOG_OGs"]

def enog_fetcher(enog_list,taxlev):
        l=[]
        for i in enog_list:
                l.append(taxlev in i)
        try:
                ans=enog_list[np.where(l)[0][0]]
        except IndexError:
                ans=None
        if ans:
                ans=ans.split("@"+taxlev+"|")[0]
        return(ans)

out_arr=np.array([])
for i in query_genes_enogs:
        out_arr=np.append(out_arr,enog_fetcher(i.split(","),"33154"))
        
enogs=np.unique(out_arr[out_arr != None])

def gene_IDs_to_enog(annot_file,gene_ID_file):
        emap_table=pd.read_csv(annot_file,sep="\t")
        if type(gene_ID_file) == str:
                inferred_ncgenes=np.genfromtxt(gene_ID_file,dtype=str)
        in_common=np.intersect1d(inferred_ncgenes,emap_table.loc[:,"#query"])
        emap_table.index=emap_table.loc[:,"#query"]
        query_genes_enogs=np.array(emap_table.loc[in_common,"eggNOG_OGs"])
        out_arr=np.array([])
        for i in query_genes_enogs:
                out_arr=np.append(out_arr,enog_fetcher(i.split(","),"33154"))
        enogs=np.unique(out_arr[out_arr != None])
        return(enogs)

import pandas as pd
import numpy as np

sp_data=pd.read_csv("Spis_adult_broad_cell_type_gene_FC.tsv",sep="\t")
nc_cells=np.genfromtxt("Sp_ncCells.list",dtype=str)
out=[]
for i in sp_data.columns:
        if i in nc_cells:
                out.append(True)
        else:
                out.append(False)
out=np.array(out)
nc_cols=np.where(out)[0]
sp_nconly=sp_data.iloc[:,nc_cols]
sp_nconly.index=sp_data.iloc[:,0] #add gene ids as index names
sp_nconly.to_csv("Sp_ncOnly_UMIcounts.csv",sep=",",columns=True,indes=True) #save

import numpy as np
import pandas as pd
#read in the UMI counts table of the neurons only
sp_ncs=pd.read_csv("Sp_ncOnly_UMIcounts.csv",sep=",",index_col=0,header=0)
#booleanize by turning to '0' all cells with a zero, and to '1' all cells with a postive count
sp_nc_bool=(sp_ncs > 0).astype(int)
# sum rows - because we're using the booleanized version of the table, this will equal to the number of cells in which each gene is found
sp_bool_rowsums=sp_nc_bool.sum(axis=1)
#get IDs for each gene at all expressed in these cells.
sp_nc_geneIDs=sp_nc_bool.index[np.where(sp_bool_rowsums > 0)[0]]
#save this list to disk
np.savetxt("Sp_genesON_Neurons.list",sp_nc_geneIDs,fmt='%s',delimiter=",")

import glob as glob
minibatch_filelist = glob.glob("train*.csv")

master=pd.DataFrame(index=None,columns=range(269),dtype=float)
col_heads=np.genfromtxt("header.tmp",delimiter=",",dtype=str)
master.columns=col_heads

in_tab=pd.read_csv("train_subset_1.csv",sep=",")
out_df = pd.concat((master,in_tab.sample(n=248832)))

del in_tab; gc.collect()

def minibatch_maker(files_list: list, chunk_size: int):
        master=pd.DataFrame(index=None,columns=range(269),dtype=float)
        col_heads=np.genfromtxt("header.tmp",delimiter=",",dtype=str)
        master.columns=col_heads
        for i in files_list:
                in_tab = pd.read_csv(i,sep=",")
                master = pd.concat((master,in_tab.sample(n=chunk_size)))
        return master

import pandas as pd
import numpy as np
# Read in MetaCell table
cteno_mcdat = pd.read_csv("Mnemiopsis_metacell_raw_counts",sep="\t",header=0)
# Make a list with the Neuronal MetaCells identified by Sachkova et al.
nc_metacells = ['27','28','29','30','31','32','33','34','35','36','40','55']
# Subsample the MetaCell table
cteno_nerv_metas = cteno_mcdat.loc[:,nc_metacells]
# extract the row indices of the genes which were expressed in any of the nerve cells at any level (extra lax)
idcs_expressed = np.where(np.array(cteno_nerv_metas.sum(axis=1) > 0))[0]
# Create the list of genes expressed at any level in any of the neuronal metacells
genes_on_neurons = np.array(cteno_nerv_metas.index[idcs_expressed])
# Save list to file
np.savetxt("Ml_inferredNeurGenes.csv", genes_on_neurons,fmt='%s',delimiter=",")


import numpy as np
import pandas as pd
import glob as glob

total_gos_fnames = glob.glob('*totalOGs.list')
nc_gos_fnames = glob.glob('*inferred_ncENOGs.list')
claves=[]
contenidos=[]
for i in total_gos_fnames:
        claves.append(i.split("_")[0])
        contenidos.append(list(np.genfromtxt(i,dtype=str)))

tot_ogs_dict = dict(zip(claves,contenidos))

import psutil
def minibatch_maker(files_list: list, chunk_size: int):
        master=pd.DataFrame(index=None,columns=range(269),dtype=float)
        col_heads=np.genfromtxt("header.tmp",delimiter=",",dtype=str)
        master.columns=col_heads
        count=1
        for i in files_list:
                in_tab = pd.read_csv(i,sep=",")
                master = pd.concat((master,in_tab.sample(n=chunk_size)))
                print(f"Subsampled {count} files - {psutil.virtual_memory().percent} of RAM used.")
                count+=1
        return master


import psutil
import numpy as np
import pandas as pd
import os

inlist=np.array([x for x in os.listdir() if "train_subset" in x])
np.random.shuffle(inlist)

import psutil
def minibatch_maker(files_list: list):
        master=pd.DataFrame(index=None,columns=range(269),dtype=float)
        col_heads=np.genfromtxt("header.tmp",delimiter=",",dtype=str)
        master.columns=col_heads
        count=1
        for i in files_list:
                in_tab = pd.read_csv(i,sep=",")
                master = pd.concat((master,in_tab))
                print(f"Subsampled {count} files - {psutil.virtual_memory().percent} of RAM used.")
                count+=1
        return master

print(f"Total amount of input files is {len(inlist)}")
v = minibatch_maker(inlist,chunk_size=10000)

import tensorflow as tf
minibatch = "2" #Something something - has to be either a pandas DF or a numpy array
minibatch = tf.data.Dataset.from_tensor_slices(minibatch)

import math
import scipy as sc
def umi_counts_transformer(in_data: pd.DataFrame):
        # Assume a UMI counts table as a numpy array.
        # check if entire dataset is of integers (even if they're formatted as float)
        if np.any(in_data.dtypes != int):
                in_data = in_data.astype(int)
        informative_rows = []
        for i in range(in_data.shape[0]):
                nonzero_idcs = in_data.iloc[i].to_numpy().nonzero()
                nonzero_vals = in_data.iloc[i].iloc[nonzero_idcs]
                # Test if counts are uninformative - if the distribution of counts resembles a uniform distribution, they're not informative.
                dummy_uniform_dist = np.repeat(np.median(nonzero_vals),len(nonzero_vals))
                if len(nonzero_vals) == 0 or sc.stats.kstest(nonzero_vals, dummy_uniform_dist).pvalue > 0.05:
                        print(f"Row {i} is not informative...skipping it")
                else:
                        print(f"Row {i} *is* informative!, recoding values.")
                        informative_rows.append(i)
                        percentile_5 = np.percentile(nonzero_vals,5)
                        percentile_95 = np.percentile(nonzero_vals,95)
                        lower_5 = np.where(in_data.iloc[i] < percentile_5)[0]
                        # Set lower 5% to 1
                        in_data.iloc[i,lower_5] = 1
                        upper_5 = np.where(in_data.iloc[i] >= percentile_95)[0]
                        # Set top 5% to 3
                        in_data.iloc[i,upper_5] = 3
                        extremes = np.hstack((lower_5,upper_5))
                        normies = np.setdiff1d(nonzero_idcs,extremes)
                        # Set all other non-zero values to 2
                        in_data.iloc[i,normies] = 2
                out_data=in_data.iloc[informative_rows]
        return out_data

def subpop_cell_selector(transformed_dataset: pd.DataFrame, markers: np.ndarray, percentile = 95):
        transf = transformed_dataset
        # Check which of the subpop markers are in the 'informative markers' dataset - exit if none
        markers_present = np.intersect1d(transf.index, markers)
        if len(markers_present) == 0:
                print("No markers are informative in this dataset - :(.")
                exit
        else:
                # Get the row index of each of the subpopulation markers, and sort them
                subpop_markers=list(map((lambda x: np.where(x == transf.index)[0][0]),markers_present))
                subpop_markers.sort()
                # transform the dataset into a boolean data frame where 'True' is values which were in 1 (downreg) OR 3 (upreg), and otherwise False
                bool_dset = np.bitwise_or((transf == 1),(transf == 3))
                # calculate the number of times markers appear on each cell
                subpop_marker_count = np.sum(bool_dset.iloc[:,subpop_markers],axis=0)
                # get the column index for all cells which are on the 95th percentile of expression of markers
                subpop_indices = np.where(subpop_marker_count >= np.percentile(subpop_marker_count,percentile))[0]
                #...add permutation test to statistically determine what other genes are preferentially expressed in the chosen subpopulation
        return(bool_dset,subpop_indices)


def diff_test(bool_data: pd.DataFrame, subpop_indices, perms = 100):
        # calculate the total of columns in the full dataset
        tot_cols = bool_data.shape[1]
        # identify the indices of the 'background' columns (those not in the subpop of interest)
        background_indices = np.setdiff1d(range(tot_cols),subpop_indices)
        # count the total of subpopulation cells for calculating proportions later
        tot_subpop = subpop_indices.size
        # count the total of background cells for calculating proportions later
        tot_background = background_indices.size
        # initialize the lists where some of the data will end up in
        percentiles = []
        obs_values = []
        for i in range(bool_data.shape[0]):
                # set up random selections of cells for permutation test - defaults to 100 permutations
                subpop_sample = np.random.choice(range(tot_cols),size=(perms,tot_subpop),replace=False) # 'subpopulation permutes'
                backg_sample = np.array(list(map((lambda x: np.setdiff1d(np.arange(tot_cols), x)), subpop_sample))) # 'background permutes'
                coords = tuple(map((lambda x,y: (x,y)),subpop_sample, backg_sample)) # tuple where [0] has the random selection that will be the 'subpopulation', and [1] has the background.
                # append the *actual* ratio of subpopulation expression for gene i
                obs_values.append((sum(bool_data.iloc[i,subpop_indices])/tot_subpop)/(sum(bool_data.iloc[i,background_indices])/tot_background))
                # Calculate the 5th and 95th percentiles for each permutation
                for_ptest = []
                for j in range(len(coords)):
                        val = expression_ratios(bool_data.iloc[i,:],coords[j])
                        for_ptest.append(val)
                percentiles.append((np.percentile(for_ptest,5),np.percentile(for_ptest,95)))
        l = list(zip(obs_values,percentiles))
        dexp_idcs = np.where(list(map((lambda x: np.bitwise_or(x[0] < x[1][0], x[0] > x[1][1])),l)))[0]
        dexp_ids = bool_data.index[dexp_idcs]
        return(dexp_ids)
        # Now, check whether obs_value[i] is <= to the percentiles[i,0] OR >= to percentiles[i,1]. BTW begin by testing the loops from above.
                
                        
                        
subpop_sample = np.random.choice(range(tot_cols),size=(100,tot_subpop),replace=False)
backg_sample = np.array(list(map((lambda x: np.setdiff1d(np.arange(tot_cols), x)), subpop_sample)))
coords = tuple(map((lambda x,y: (x,y)),subpop_sample, backg_sample))

def expression_ratios(series,coords):
        qty1 = len(coords[0])
        qty2 = len(coords[1])
        return (sum(series.iloc[coords[0]]/qty1))/(sum(series.iloc[coords[1]]/qty2))

percentiles = []
for i in range(bool_data.shape[0]):
        l = []
        for j in range(len(coords)):
               val = expression_ratios(bool_data.iloc[i,:],coords[j])
               print(f"Row {i}, permutation {j}: {val}") # DEBUG
               l.append(val)
        l = np.array(l)
        percentiles.append((np.percentile(l,5),np.percentile(l,95)))



(lambda x: np.bitwise_or(x[0] < x[1][0], x[0] > x[1][1]))
dexp_idcs = np.where(list(map((lambda x: np.bitwise_or(x[0] < x[1][0], x[0] > x[1][1])),l)))[0]
dexp_ids = bool_data.index[dexp_idcs]

import numpy as np
import pandas as pd
import scanpy as sp
import scipy as sc
dm_ad = sp.read_loom("Dm_scRNAseq.loom")
dm_data = sc.sparse.csr_matrix.todense(dm_ad.X)
dm_data = dm_data.T

dm_df = pd.DataFrame(data=dm_data, index=dm_ad.var.Accession, columns=dm_ad.obs.index)
markers = np.genfromtxt("NeuralMarkers_fly.csv",dtype=str)


import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sp
import scipy as sc

dr_ad=sp.read_h5ad("Dr_scRNAseq_counts.h5ad")
dr_data = dr_ad.X.todense()
dr_data = dr_data.astype(int) # Took longer than expected
dr_data = dr_data.T # Took shorter than expected
rownames=dr_ad.var.gene_ids
colnames = dr_ad.obs.index
dr_pddata = pd.DataFrame(dr_data, index = rownames, columns = colnames, dtype = int)
dr_pddata.to_csv("Dr_scRNAseq-cts.csv",sep=",",header=True, index = True) # Will this ever be done?


def minibatch_maker(files_list):
        count=0
        for i in np.arange(len(files_list)):
                curr_file = files_list[i]
                if i == 0:
                        master = np.genfromtxt(curr_file, delimiter=",", skip_header=0)
                else:
                        master = np.concatenate((master,np.genfromtxt(curr_file, delimiter=",", skip_header=0)))
                        print(f"Subsampled {count+1} files - {psutil.virtual_memory().percent}% of RAM used.")
        return(master)


# CHATGPT CODE TO BE MODIFIED for CtenoNN training
import tensorflow as tf
import numpy as np

class BigTrainingDataLoader:
    def __init__(self, data, batch_size):
        self.data = data
        #self.labels = labels # I think I don't need this.
        self.batch_size = batch_size
        self.num_samples = len(data) # list of snippett filenames
        self.num_batches = int(np.ceil(self.num_samples / self.batch_size))
        self.current_batch = 0

    def next_batch(self):
        # Batch Size gets set upon running .__init__, and must be the number of full datasets which fit into the RAM.
        # In this function I have to set up the loop that reads adds into the batch as many snippetts as possible, and outputs the data and labels rows separately (see blow)
        # Assuming I know what the 'batch_size' is, I can already start!
        start_idx = self.current_batch * self.batch_size
        end_idx = min((self.current_batch + 1) * self.batch_size, self.num_samples)

        batch_data = self.data[start_idx:end_idx]
        #batch_labels = self.labels[start_idx:end_idx] # I load data and labels everything as one dataset.

        # Update the current batch index for the next iteration
        self.current_batch = (self.current_batch + 1) % self.num_batches
        # output here has to separate the data (cols 0:ncol(data)-1 for 'batch_data') from the labels (last col - 'batch_labels')
        return batch_data, batch_labels

    def reset(self):
        self.current_batch = 0

# Example usage:
# Assume you have your data and labels as numpy arrays
data = np.random.randn(100, 32)  # Example data with 100 samples of size 32
labels = np.random.randint(0, 2, size=(100,))  # Example labels (binary classification)

batch_size = 32
data_loader = BigTrainingDataLoader(data, labels, batch_size)

# Iterate through epochs and batches
num_epochs = 3
for epoch in range(num_epochs):
    data_loader.reset()  # Reset the batch index at the beginning of each epoch

    for _ in range(data_loader.num_batches):
        batch_data, batch_labels = data_loader.next_batch()

        # Perform your training or evaluation with the current batch
        # For example, you can feed this batch to your model and update the parameters

        # Replace the following with your actual training code
        # model.train_on_batch(batch_data, batch_labels)
        

import scanpy as sp
import pandas as pd
import numpy as np

adata = sp.read_loom("Dm_scRNAseq.loom")
ns_cells = np.where(np.bitwise_or(adata.obs['R_annotation_broad'] == "neuron",adata.obs['R_annotation_broad'] == "sensory neuron"))[0]
ns_cell_ids = adata.obs.index[ns_cells]
np.savetxt("Dm_ncCells.csv",ns_cell_ids,delimiter=",",fmt='%s')

cns_cells = np.where(adata.obs['zebrafish_anatomy_ontology_class'] == 'central_nervous_system')[0]
cns_cell_bcode = adata.obs.index[cns_cells]

adata = sp.read_h5ad("Dr_scRNAseq.h5ad")
cns_cells = np.where(adata.obs['zebrafish_anatomy_ontology_class'] == 'central_nervous_system')[0]
cns_cell_bcode = adata.obs.index[cns_cells]
transf = pd.read_csv("fish_transf_table.csv",header=0, index_col=0, nrows = 100)

out_list = []
for i in cns_cell_bcode:
        val=np.where(i == transf.columns)[0]
        print(f"{i} matches to column {val}")
        out_list.append(val)
        
subpop_indices = np.hstack(out_list)
np.savetxt("fish_subpop_indices.txt",subpop_indices,dtype="%s")

import percentile_sorter as ps
np = ps.np
pd = ps.pd
org = "fly"

indices_file = org+"_subpop_indices.txt"
subpop_indices = np.genfromtxt(indices_file, dtype=int)

transf_file = org+"_transf_table.csv"
transf = pd.read_csv(transf_file, header=0, index_col = 0)
bool_dset = np.bitwise_or((transf == 1),(transf == 3))
de_gids = ps.diff_test(bool_dset, subpop_indices, 100, 5)

np.savetxt(org+"_NSGenes.list",de_gids, fmt='%s')

def percent_shared(key1, key2):
        set1=d[key1]
        set2=d[key2]
        sizes = (set1.size, set2.size)
        max_one,min_one = (set1, set2)[np.argmax(sizes)],(set1,set2)[np.setdiff1d((0,1),np.argmax(sizes))[0]]
        shared_vals=np.intersect1d(max_one,min_one)
        perc_shared_max=100*shared_vals.size/max_one.size
        perc_shared_min=100*shared_vals.size/min_one.size
        print(f"{key1}: {max(sizes)} nsENOGs\n{key2}: {min(sizes)} nsENOGs\n {perc_shared_max}, and {perc_shared_min}")
        return
#################################################

def one_permutation(ns_a, ns_b, tot_a,tot_b):
        rng = np.random.default_rng()
        perm_a = rng.choice(tot_a,size=ns_a.size, replace=False)
        perm_b = rng.choice(tot_b,size=ns_b.size, replace=False)
        return(compare_shared(perm_a,perm_b))

def compare_shared(set1,set2):
        shared=np.round(np.intersect1d(set1,set2).size,1)
        return(shared)

#same as above, but assumes a sample size of 100 and only returns the shared portion:
def compare_shared(set1,set2):
        shared=np.intersect1d(set1,set2).size
        return(shared)

ssize=ns1.size

def characterize_vector(vect):
        print(f"Min: {vect.min()}\nMedian: {np.median(vect)}\nMax: {vect.max()}")

np.array([ one_permutation(ns1, ns2, tot1, tot2) for i in range(100) ])


####### FORMAL SCRIPT
import numpy as np
import matplotlib.pyplot as plt

def compare_shared(set1,set2,curr_scale):
        shared=np.round(np.intersect1d(set1,set2).size*curr_scale*100,1)
        return(shared)

def one_permutation(ns_a, ns_b, tot_a,tot_b, curr_scale):
        rng = np.random.default_rng()
        perm_a = rng.choice(tot_a,size=ns_a.size, replace=False)
        perm_b = rng.choice(tot_b,size=ns_b.size, replace=False)
        return(compare_shared(perm_a,perm_b, curr_scale))

scale_count = 1826 # Nematostella's count of nsENOGs (the one with the least nsENOGs)

sp1_file = "sp1.csv" # get the filename of the first file from the args, or standardize it to 'sp1.csv' and 'spp2.csv'
sp2_file = "spp2.csv"
name1 = np.genfromtxt(sp1_file, dtype=str, delimiter=",")
names2 = np.genfromtxt(sp2_file, dtype=str, delimiter=",")
prefix1 = name1.split(" ")[0][0]+name1.split(" ")[1][0] # this will prefix the sp1 filenames.
ns1 = np.genfromtxt(prefix1+"_nsENOGs.list",dtype=str)
tot1 = np.genfromtxt(prefix1+"_totENOGs.list",dtype=str)
for i in range(len(names2)):
        name2 = names2[i]
        prefix2 = name2.split("_")[0][0]+name2.split("_")[1][0] # this will prefix the sp2 filenames, note I've used 'i' so I can then transform it into a loop easily.
        ns2 = np.genfromtxt(prefix2+"_nsENOGs.list",dtype=str)
        curr_scale = scale_count/ns2.size
        tot2 = np.genfromtxt(prefix2+"_totENOGs.list",dtype=str)
        obs = compare_shared(ns1, ns2, curr_scale)
        perms_vals = np.array([ one_permutation(ns1, ns2, tot1, tot2, curr_scale) for i in range(100) ])
        outfile_name = prefix1+"_vs_"+prefix2+"_hist.png"
        ax = plt.gca()
        ax.set_xlim([0,100])
        plt.hist(perms_vals, bins=10)
        plt.axvline(obs, color="r")
        plt.title("Expected % of shared OGs between "+ str(name1) + " and " + str(name2) + " (blue, N=100)\nwith observed % of shared nervous system OGs (red)", loc='left', fontsize='small')
        plt.xlabel("% shared OGs")
        plt.ylabel("Counts")
        plt.savefig(outfile_name)
        plt.close()

import glob
import numpy as np
import pandas as pd

tot_files=glob.glob("*totENOGs.list")
ref_file="eumetazoa_nsENOGs.list"

d = {}

for i in tot_files:
        key = i.split("_")[0]
        val = np.genfromtxt(i, dtype=str)
        d[key] = val
        
ref = np.genfromtxt(ref_file, dtype=str)

def compare_shared(set1,set2,curr_scale):
        #### WATCH OUT!!! % calculated over the second set's total size.
        shared=np.round(np.intersect1d(set1,set2).size/set2.size*curr_scale*100,1)
        return(shared)
ids = []
perc_shared = []

for i in d:
        ids.append(i)
        perc_shared.append(compare_shared(d[i], ref))

# Now we put that all into a nice pandas DF...
tot_comparison=pd.DataFrame(perc_shared, index=ids)
tot_comparison.columns = ["% eumetazoa nsENOGs in genome"]
# and we reorder the indices IDs to more or less follow a phylogenetic branching pattern.
reorder = ["Sr","Mb","Sl","Em","Aq","Ml","Pb","Ta","Nv","Hv","Sp","Sm","Sa","Cg","Dm","Su","Cc","Dr", "Xt","Mm"]
cmp_table = tot_comparison.reindex(reorder)
# Should I add also a 'N/M' string? Nah, that would give too much focus on the specific ENOGs, and I'd prefer to avoid that.
# Time to save it!
cmp_table.to_csv("EumetNSENOGs_genome_representation.csv",sep=",")

# AI-generated code to infer a phylogenetic tree using ML
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo import PhyloML

def estimate_phylogenetic_tree(alignment_file):
    # Load the alignment from a file
    alignment = AlignIO.read(alignment_file, "fasta")

    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Construct the tree using the maximum likelihood method
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)

    # Perform maximum likelihood optimization
    tree = PhyloML(tree, alignment) # check if it explores different tree topologies - otherwise consider several independent random starting trees or smth

    return tree

# Example usage
alignment_file = "example_alignment.fasta"
phylogenetic_tree = estimate_phylogenetic_tree(alignment_file)

# Print the phylogenetic tree
Phylo.draw(phylogenetic_tree)

import numpy as np
import pandas as pd
# First, build a dictionary by importing each of the species' prefixes as the key, and the nsENOG list as the value, let's call it 'd'
d ={'Cg': np.genfromtxt("Crassostrea_gigas/Cg_nsENOGs.list",dtype=str), 'Dr': np.genfromtxt("Danio_rerio/Dr_inf_nsENOGs.list", dtype=str),
    'Dm': np.genfromtxt("Drosophila_melanogaster/Dm_inf_nsENOGs.list", dtype=str),
    'Ml': np.genfromtxt("Mnemiopsis_leidyi/Ml_nsENOGs.list",dtype=str),
    'Mm': np.genfromtxt("Mus_musculus/Mm_inf_nsENOGs.list",dtype=str) }
# Automate this, don't type it!

l=[] # a list of list with all the nsENOGs
for i in d:
        l.append(list(d[i]))

# flatten and deduplicate the list
ddl=np.unique([item for sublist in l for item in sublist])

def where_is_enog(enog_id, dix):
        l=[]
        for i in dix:
                l.append(int(enog_id in dix[i]))
        return(l)

# get the results
x = np.array([ where_is_enog(i, d) for i in ddl ])
data_table=pd.DataFrame(x, index=ddl, columns=d.keys()) # make the table
data_table.sum(axis=1) # add numbers by row, the columns in which the sum is = to the number of columns, should be 39

def shared_enogs(species: list, data_table: pd.DataFrame):
        tot = len(species)
        shared_enogs=np.array(data_table.index[data_table[species].sum(axis=1) == tot])
        return(shared_enogs)

# make a list of all the possible species combinations
from itertools import combinations
keys_set=set(d.keys())
combi_list=[]
for n in range(2,len(keys_set)+1):
        curr_comb=list(combinations(keys_set,n))
        for i in curr_comb:
                combi_list.append(list(i))

# combi_list has all the possible combinations which can be given to the 'shared_enogs' function.      
# that loop produces the counts of shared ENOGs
# Now, let's format these results into a pd.DataFrame for a bar graph
# First, the labels:
idcs=[ "_".join(i) for i in combi_list ]
# Now the values:
vals=list(map((lambda x: shared_enogs(x, data_table).size),combi_list))
combis_enog_counts=pd.DataFrame(vals, index=idcs)
combis_enog_counts.columns=["shared ENOGs"]
reordered=combis_enog_counts.sort_values("shared ENOGs",ascending=False)

ax = reordered.plot(kind="bar", title='Shared nsENOGs per species combinations', xlabel='Organism set', xticks=[], ylabel='count shared ENOGs', width=width, figsize=(8,8))
width=0.5
labels=[ x.replace("_",", ") for x in reordered.index]
###
def autolabel(labels,width):
        for i in range(len(labels)):
                ax.text(i-(width/2), reordered.iloc[i,0], labels[i], ha='left', va='bottom', rotation=80., fontsize='small')
autolabel(labels, width)
plt.savefig("Bar_graph.png")
plt.close()


import pandas as pd

table = pd.read_csv("Out_Binary_table.csv", sep=",", header=0, index_col=0)

for_fasta=table.iloc[:,0:7].T
for_fasta.to_csv("nsENOGs_presabs.fas", sep=",")